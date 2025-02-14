import os
import sys
from pathlib import Path

# Add the project root directory to Python path
project_root = str(Path(__file__).parent.parent.parent)
if project_root not in sys.path:
    sys.path.append(project_root)

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import seaborn as sns
import logging
import numpy as np
from typing import Dict, List, Tuple, Optional
from cyvcf2 import VCF
import pandas as pd
import subprocess
from datetime import datetime
from dask import delayed, compute

from src.utils.hpc import HPCManager

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class GenomeScorer:
    """
    Enhanced version of SprintGeneticsAnalyzer with additional scoring capabilities
    and HPC optimization using Dask.
    """
    
    def __init__(self, config_path: Optional[Path] = None):
        self.snp_data = None
        self.scores = {}
        self.vcf_path = None
        self.config = self._load_config(config_path) if config_path else {}
        self.hpc_manager = None
        self.initialization_status = {
            'config_loaded': False,
            'hpc_initialized': False,
            'data_loaded': False
        }
        
        try:
            self.hpc_manager = HPCManager()
            self.hpc_manager.start_cluster()
            self.initialization_status['hpc_initialized'] = True
        except Exception as e:
            logger.warning(f"HPC initialization failed: {e}. Will use sequential processing.")
        
        # Updated sprint variants based on your genome
        self.sprint_variants = {
            'INSIG2': {
                'rs7566605': {
                    'effect': 'Body composition and muscle development',
                    'beneficial_allele': 'G',
                    'risk_allele': 'C',
                    'impact': 'Influences body mass index and muscle composition'
                }
            },
            'ACTN3': {
                'rs1815739': {
                    'effect': 'Speed/power performance',
                    'beneficial_allele': 'C',
                    'risk_allele': 'T',
                    'impact': 'R577X variant affects fast-twitch muscle fiber production'
                }
            },
            'BDKRB2': {
                'rs1799722': {
                    'effect': 'Exercise efficiency',
                    'beneficial_allele': 'C',
                    'risk_allele': 'T',
                    'impact': 'Influences exercise efficiency and performance'
                }
            },
            'EPAS1': {
                'rs1867785': {
                    'effect': 'Oxygen processing',
                    'beneficial_allele': 'G',
                    'risk_allele': 'A',
                    'impact': 'Affects hemoglobin concentration and oxygen processing'
                }
            },
            'GNPDA2': {
                'rs10938397': {
                    'effect': 'Energy metabolism',
                    'beneficial_allele': 'A',
                    'risk_allele': 'G',
                    'impact': 'Influences body composition and energy utilization'
                }
            },
            'CCL2': {
                'rs2857656': {
                    'effect': 'Recovery and inflammation',
                    'beneficial_allele': 'G',
                    'risk_allele': 'C',
                    'impact': 'Affects post-exercise inflammation and recovery'
                }
            },
            'NOS3': {
                'rs2070744': {
                    'effect': 'Blood flow and oxygen delivery',
                    'beneficial_allele': 'T',
                    'risk_allele': 'C',
                    'impact': 'Influences nitric oxide production and blood flow'
                }
            },
            'CKM': {
                'rs8111989': {
                    'effect': 'Energy production in muscles',
                    'beneficial_allele': 'C',
                    'risk_allele': 'T',
                    'impact': 'Affects creatine kinase activity and energy availability'
                }
            },
            'IL6': {
                'rs1800795': {
                    'effect': 'Recovery and inflammation',
                    'beneficial_allele': 'G',
                    'risk_allele': 'C',
                    'impact': 'Affects post-exercise recovery and adaptation'
                }
            },
            'MSTN': {
                'rs1805086': {
                    'effect': 'Muscle growth and development',
                    'beneficial_allele': 'C',
                    'risk_allele': 'T',
                    'impact': 'Influences muscle mass and strength development'
                }
            },
            'PPARGC1A': {
                'rs8192678': {
                    'effect': 'Energy metabolism',
                    'beneficial_allele': 'G',
                    'risk_allele': 'A',
                    'impact': 'Influences mitochondrial function and energy production'
                }
            },
            'UCP2': {
                'rs659366': {
                    'effect': 'Energy efficiency',
                    'beneficial_allele': 'T',
                    'risk_allele': 'C',
                    'impact': 'Influences metabolic efficiency during exercise'
                }
            },
            'COMT': {
                'rs4680': {
                    'effect': 'Pain tolerance and recovery',
                    'beneficial_allele': 'G',
                    'risk_allele': 'A',
                    'impact': 'Affects pain perception and recovery capacity'
                }
            },
            'CHRM2': {
                'rs324640': {
                    'effect': 'Cardiac function',
                    'beneficial_allele': 'G',
                    'risk_allele': 'A',
                    'impact': 'Influences heart rate response to exercise'
                }
            }
        }
        
    def _load_config(self, config_path: Path) -> Dict:
        """Load configuration for scoring weights and parameters."""
        try:
            config = pd.read_json(config_path).to_dict()
            # Merge with default sprint variants if provided in config
            if 'sprint_variants' in config:
                self.sprint_variants.update(config['sprint_variants'])
            return config
        except Exception as e:
            logger.error(f"Error loading config: {e}")
            return {}
            
    def get_status(self) -> Dict:
        """Get current status of the scorer."""
        return {
            'initialization': self.initialization_status,
            'variants_loaded': len(self.snp_data) if self.snp_data is not None else 0,
            'scores_computed': len(self.scores) > 0,
            'hpc_available': self.hpc_manager is not None and self.hpc_manager.client is not None
        }

    def load_data(self, vcf_path: Path) -> Dict:
        """Load SNP data from VCF file with detailed logging."""
        self.vcf_path = vcf_path
        status = {'success': False, 'variants_found': 0, 'errors': [], 'details': []}
        
        try:
            if not (str(vcf_path).endswith('.vcf.gz') or str(vcf_path).endswith('.snp.vcf.gz')):
                logger.error("VCF file must be gzipped (.vcf.gz or .snp.vcf.gz)")
                raise ValueError("VCF file must be gzipped (.vcf.gz or .snp.vcf.gz)")

            # Create index if needed
            self._create_tabix_index(vcf_path)

            # Now load the VCF
            logger.info(f"Opening VCF file: {vcf_path}")
            vcf_reader = VCF(str(vcf_path))
            variants_data = []
            
            # Log some basic VCF info
            logger.info(f"VCF samples: {vcf_reader.samples}")
            logger.info(f"Number of variants to search for: {sum(len(v) for v in self.sprint_variants.values())}")
            
            for gene, variants in self.sprint_variants.items():
                for rsid, info in variants.items():
                    try:
                        logger.info(f"Searching for variant {rsid} (gene: {gene})")
                        found = False
                        
                        # Try different ways to query the variant
                        # Method 1: Direct rsID
                        try:
                            for record in vcf_reader(rsid):
                                genotype = self._get_genotype(record)
                                variants_data.append({
                                    'variant': rsid,
                                    'gene': gene,
                                    'genotype': genotype,
                                    'chromosome': record.CHROM,
                                    'position': record.POS
                                })
                                found = True
                                logger.info(f"Found {rsid} - Genotype: {genotype}")
                                break
                        except Exception as e:
                            logger.debug(f"Direct rsID lookup failed for {rsid}: {e}")
                        
                        # Method 2: Try without 'rs' prefix
                        if not found and rsid.startswith('rs'):
                            try:
                                variant_id = rsid[2:]  # Remove 'rs' prefix
                                for record in vcf_reader(variant_id):
                                    genotype = self._get_genotype(record)
                                    variants_data.append({
                                        'variant': rsid,
                                        'gene': gene,
                                        'genotype': genotype,
                                        'chromosome': record.CHROM,
                                        'position': record.POS
                                    })
                                    found = True
                                    logger.info(f"Found {rsid} (without rs prefix) - Genotype: {genotype}")
                                    break
                            except Exception as e:
                                logger.debug(f"Non-rs lookup failed for {rsid}: {e}")
                        
                        if not found:
                            status['details'].append(f"Variant {rsid} not found in VCF")
                            logger.warning(f"Could not find variant {rsid} in VCF file")
                            
                    except Exception as e:
                        error_msg = f"Error processing {rsid}: {str(e)}"
                        status['errors'].append(error_msg)
                        logger.error(error_msg)
                
            self.snp_data = pd.DataFrame(variants_data)
            status['variants_found'] = len(self.snp_data)
            status['success'] = True
            self.initialization_status['data_loaded'] = True
            
            # Log summary
            logger.info(f"Found {status['variants_found']} variants out of {sum(len(v) for v in self.sprint_variants.values())} searched")
            if status['variants_found'] > 0:
                logger.info("Found variants:")
                for _, row in self.snp_data.iterrows():
                    logger.info(f"  {row['gene']} - {row['variant']}: {row['genotype']}")
            
        except Exception as e:
            error_msg = f"Critical error loading VCF: {str(e)}"
            status['errors'].append(error_msg)
            logger.error(error_msg)
            
        return status

    def _create_tabix_index(self, vcf_path: Path) -> None:
        """Create tabix index for VCF file if it doesn't exist."""
        try:
            # Check if index exists
            if not os.path.exists(str(vcf_path) + '.tbi'):
                logger.info(f"Creating tabix index for {vcf_path}")
                
                # Create index using tabix
                subprocess.run(['tabix', '-p', 'vcf', str(vcf_path)], check=True)
                logger.info("Tabix index created successfully")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error creating tabix index: {e}")
            raise
        except FileNotFoundError:
            logger.error("tabix command not found. Please install tabix (part of htslib)")
            raise

    def _get_genotype(self, record) -> str:
        """Extract genotype from VCF record with better error handling."""
        try:
            # Get the genotype for the first sample (assuming single-sample VCF)
            gt = record.genotypes[0]
            
            # Convert numeric genotype to alleles
            alleles = [record.ALT[i-1] if i > 0 else record.REF for i in gt[:2]]
            genotype = ''.join(alleles)
            
            logger.debug(f"Raw genotype: {gt}, Converted to: {genotype}")
            return genotype
            
        except Exception as e:
            logger.error(f"Error extracting genotype: {e}")
            return "NN"  # Return placeholder for failed genotype extraction

    def _analyze_genotype(self, genotype: str, beneficial: str, risk: str) -> bool:
        """Analyze if genotype is beneficial for sprint performance."""
        if not genotype:
            return False
        
        alleles = genotype.replace('|', '/').split('/')
        beneficial_count = sum(1 for allele in alleles if allele == beneficial)
        return beneficial_count >= 1

    @delayed
    def _score_variant(self, variant: str, genotype: str) -> float:
        """Score individual variant (can be executed in parallel)."""
        try:
            for gene, variants in self.sprint_variants.items():
                if variant in variants:
                    info = variants[variant]
                    return self._calculate_variant_score(genotype, 
                                                      info['beneficial_allele'], 
                                                      info['risk_allele'])
        except Exception as e:
            logger.warning(f"Error scoring variant {variant}: {e}")
            return 0.0
        return 0.0

    def calculate_sprint_score(self) -> Dict:
        """Calculate sprint performance score with detailed status reporting."""
        if self.snp_data is None:
            return {
                'success': False,
                'error': 'No data loaded',
                'partial_results': {}
            }

        results = {
            'success': False,
            'total_score': 0.0,
            'variant_scores': {},
            'errors': [],
            'processed_variants': 0
        }

        try:
            computed_scores = []
            
            if self.hpc_manager and self.hpc_manager.client:
                # Parallel processing with HPC
                try:
                    with self.hpc_manager.client:
                        tasks = []
                        for _, row in self.snp_data.iterrows():
                            task = self._score_variant(row['variant'], row['genotype'])
                            tasks.append(task)
                        computed_scores = compute(*tasks)
                except Exception as e:
                    results['errors'].append(f"HPC processing failed: {e}")
                    # Fallback to sequential processing
                    logger.warning("Falling back to sequential processing")
                    computed_scores = [
                        self._score_variant(row['variant'], row['genotype']).compute()
                        for _, row in self.snp_data.iterrows()
                    ]
            else:
                # Sequential processing
                computed_scores = [
                    self._score_variant(row['variant'], row['genotype']).compute()
                    for _, row in self.snp_data.iterrows()
                ]

            # Process results
            results['variant_scores'] = dict(zip(self.snp_data['variant'], computed_scores))
            results['total_score'] = np.mean(computed_scores) if computed_scores else 0.0
            results['processed_variants'] = len(computed_scores)
            results['success'] = True

        except Exception as e:
            results['errors'].append(f"Score calculation error: {str(e)}")
            
        return results

    def generate_report(self, output_dir: Path = Path('results')) -> None:
        """Generate comprehensive report with visualizations."""
        if not self.scores:
            logger.error("No scores available. Run calculate_sprint_score first.")
            return

        output_dir.mkdir(exist_ok=True)
        
        # Generate visualizations
        self._plot_variant_distribution(output_dir)
        self._plot_performance_score(output_dir)
        self._generate_text_report(output_dir)

    def _plot_variant_distribution(self, output_dir: Path) -> None:
        """Create a bar plot of variant distribution."""
        plt.figure(figsize=(12, 6))
        
        data = {
            'Beneficial': self.scores['summary']['beneficial_variants'],
            'Risk': self.scores['summary']['risk_variants']
        }
        
        sns.barplot(x=list(data.keys()), y=list(data.values()))
        plt.title('Distribution of Sprint-Related Genetic Variants')
        plt.ylabel('Number of Variants')
        
        plt.savefig(output_dir / 'variant_distribution.png')
        plt.close()

    def _plot_performance_score(self, output_dir: Path) -> None:
        """Create a gauge plot of overall performance score."""
        fig, ax = plt.subplots(figsize=(10, 10))
        
        score = self.scores['total_score']
        
        # Create gauge
        ax.set_theta_direction(-1)
        ax.set_theta_offset(np.pi/2.0)
        
        bounds = np.array([0, 20, 40, 60, 80, 100])
        theta = np.linspace(0, 2*np.pi, len(bounds))
        
        ax.plot(theta, bounds, color='k', linewidth=2)
        ax.fill(theta, bounds, score, alpha=0.3)
        
        ax.set_title(f'Sprint Performance Score: {score:.1f}%')
        
        plt.savefig(output_dir / 'performance_score.png')
        plt.close()

    def _generate_text_report(self, output_dir: Path) -> None:
        """Generate a detailed text report."""
        report = [
            "Sprint Genetics Analysis Report",
            "============================\n",
            f"Overall Performance Score: {self.scores['total_score']:.1f}%\n",
            "Variant Analysis:",
            "---------------"
        ]
        
        for variant, score in self.scores['variant_scores'].items():
            gene_info = next((
                (gene, info) 
                for gene, variants in self.sprint_variants.items() 
                for v, info in variants.items() 
                if v == variant
            ), (None, None))
            
            if gene_info[0]:
                report.extend([
                    f"\nGene: {gene_info[0]}",
                    f"Variant: {variant}",
                    f"Effect: {gene_info[1]['effect']}",
                    f"Impact: {gene_info[1]['impact']}",
                    f"Status: {'Beneficial' if score > 0 else 'Risk'} variant"
                ])
        
        report.extend([
            "\nSummary:",
            "--------",
            f"Total variants analyzed: {self.scores['summary']['total_variants']}",
            f"Beneficial variants: {self.scores['summary']['beneficial_variants']}",
            f"Risk variants: {self.scores['summary']['risk_variants']}"
        ])
        
        with open(output_dir / 'sprint_genetics_report.txt', 'w') as f:
            f.write('\n'.join(report))

    def __del__(self):
        """Cleanup method."""
        if self.hpc_manager:
            try:
                self.hpc_manager.stop_cluster()
            except:
                pass
