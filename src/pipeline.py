import os
import sys
from pathlib import Path
import logging
from typing import Dict, Optional, Union
import json
import time
from datetime import datetime
import pandas as pd

# Add the project root directory to Python path
project_root = str(Path(__file__).parent.parent)
if project_root not in sys.path:
    sys.path.append(project_root)

from src.core.scoring import GenomeScorer
from src.core.network import NetworkBuilder
from src.core.database import DatabaseIntegrator
from src.utils.hpc import HPCManager

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class SprintGenomePipeline:
    """
    Main pipeline class that orchestrates the entire analysis workflow.
    Manages the execution of genome scoring, network analysis, and database integration.
    """
    
    def __init__(
        self,
        vcf_path: Union[str, Path],
        output_dir: Union[str, Path],
        config_path: Optional[Union[str, Path]] = None,
        use_hpc: bool = True
    ):
        self.vcf_path = Path(vcf_path)
        self.output_dir = Path(output_dir)
        self.config_path = Path(config_path) if config_path else None
        self.use_hpc = use_hpc
        
        # Create output directories
        self.results_dir = self.output_dir / 'results'
        self.cache_dir = self.output_dir / 'cache'
        self.network_dir = self.output_dir / 'networks'
        self.reports_dir = self.output_dir / 'reports'
        
        for directory in [self.results_dir, self.cache_dir, 
                         self.network_dir, self.reports_dir]:
            directory.mkdir(parents=True, exist_ok=True)
            
        # Initialize HPC manager first
        self.hpc_manager = None
        if use_hpc:
            self.hpc_manager = HPCManager()
            self.hpc_manager.start_cluster()  # Ensure cluster is started
            
        # Initialize components with HPC manager
        self.genome_scorer = GenomeScorer(config_path=self.config_path)
        self.db_integrator = DatabaseIntegrator(cache_dir=self.cache_dir)
        self.network_builder = None  # Will be initialized after scoring
        
        # Store results
        self.results = {}
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
    def run_pipeline(self) -> Dict:
        """
        Execute the complete analysis pipeline.
        Returns a dictionary containing all results.
        """
        try:
            logger.info("Starting Sprint Genome Analysis Pipeline")
            start_time = time.time()
            
            # Step 1: Load and score genetic variants
            self._run_genome_scoring()
            
            # Step 2: Integrate database information
            self._run_database_integration()
            
            # Step 3: Build and analyze network
            self._run_network_analysis()
            
            # Step 4: Generate reports
            self._generate_reports()
            
            # Save final results
            self._save_results()
            
            execution_time = time.time() - start_time
            logger.info(f"Pipeline completed successfully in {execution_time:.2f} seconds")
            
            return self.results
            
        except Exception as e:
            logger.error(f"Pipeline execution failed: {e}")
            raise
        finally:
            if self.hpc_manager and self.hpc_manager.client:
                self.hpc_manager.stop_cluster()
                
    def _run_genome_scoring(self) -> None:
        """Execute genome scoring component."""
        logger.info("Starting genome scoring analysis")
        
        try:
            # Load VCF data
            self.genome_scorer.load_data(self.vcf_path)
            
            # Calculate scores
            scores = self.genome_scorer.calculate_sprint_score()
            
            # Store results
            self.results['genome_scoring'] = scores
            
            # Generate initial reports
            self.genome_scorer.generate_report(self.reports_dir)
            
            logger.info("Genome scoring completed successfully")
            
        except Exception as e:
            logger.error(f"Genome scoring failed: {e}")
            raise
            
    def _run_database_integration(self) -> None:
        """Execute database integration component."""
        logger.info("Starting database integration")
        
        try:
            # Get genes from scoring results
            if not self.results.get('genome_scoring'):
                logger.warning("No genome scoring results available for database integration")
                return
                
            genes = list(self.results['genome_scoring'].get('variant_scores', {}).keys())
            
            if not genes:
                logger.warning("No genes available for database integration")
                return
                
            # Enrich gene data
            if self.use_hpc and self.hpc_manager and self.hpc_manager.client:
                with self.hpc_manager.client:
                    enriched_data = self.db_integrator.enrich_gene_data(genes)
            else:
                # Fallback to non-parallel processing
                logger.info("Using sequential processing for database integration")
                enriched_data = self.db_integrator.enrich_gene_data(genes)
                
            # Store results
            self.results['database_integration'] = enriched_data
            
            logger.info("Database integration completed successfully")
            
        except Exception as e:
            logger.error(f"Database integration failed: {e}")
            raise
            
    def _run_network_analysis(self) -> None:
        """Execute network analysis component."""
        logger.info("Starting network analysis")
        
        try:
            # Initialize network builder with combined data
            self.network_builder = NetworkBuilder(
                genome_data=self.results['genome_scoring'],
                config_path=self.config_path
            )
            
            # Build network
            network = self.network_builder.build_network()
            
            # Analyze network
            network_analysis = self.network_builder.analyze_network()
            
            # Generate network visualization
            self.network_builder.visualize_network(self.network_dir)
            
            # Store results
            self.results['network_analysis'] = network_analysis
            
            logger.info("Network analysis completed successfully")
            
        except Exception as e:
            logger.error(f"Network analysis failed: {e}")
            raise
            
    def _generate_reports(self) -> None:
        """Generate comprehensive reports from all analyses."""
        logger.info("Generating final reports")
        
        try:
            # Create summary report
            summary = {
                'execution_timestamp': self.timestamp,
                'genome_scoring_summary': self.results['genome_scoring']['summary'],
                'network_summary': self.results['network_analysis']['summary'],
                'database_integration_summary': {
                    'genes_analyzed': len(self.results['database_integration']['protein_interactions']),
                    'pathways_found': len(self.results['database_integration']['pathways'])
                }
            }
            
            # Save summary report
            summary_path = self.reports_dir / f'summary_report_{self.timestamp}.json'
            with open(summary_path, 'w') as f:
                json.dump(summary, f, indent=4)
                
            # Create detailed Excel report
            self._generate_excel_report()
            
            logger.info("Report generation completed successfully")
            
        except Exception as e:
            logger.error(f"Report generation failed: {e}")
            raise
            
    def _generate_excel_report(self) -> None:
        """Generate detailed Excel report with multiple sheets."""
        excel_path = self.reports_dir / f'detailed_report_{self.timestamp}.xlsx'
        
        with pd.ExcelWriter(excel_path) as writer:
            # Genome scoring results
            pd.DataFrame(self.results['genome_scoring']['variant_scores'].items(),
                        columns=['Variant', 'Score']).to_excel(writer, 
                        sheet_name='Variant_Scores', index=False)
                        
            # Network analysis results
            pd.DataFrame(self.results['network_analysis']['centrality']['degree'].items(),
                        columns=['Gene', 'Centrality']).to_excel(writer,
                        sheet_name='Network_Centrality', index=False)
                        
            # Database integration results
            pd.DataFrame(self.results['database_integration']['pathways']).to_excel(
                writer, sheet_name='Pathway_Analysis', index=False)
                
    def _save_results(self) -> None:
        """Save all results to disk."""
        logger.info("Saving pipeline results")
        
        try:
            # Save complete results as JSON
            results_path = self.results_dir / f'complete_results_{self.timestamp}.json'
            
            # Convert non-serializable objects to strings
            serializable_results = self._prepare_results_for_serialization(self.results)
            
            with open(results_path, 'w') as f:
                json.dump(serializable_results, f, indent=4)
                
            logger.info(f"Results saved to {results_path}")
            
        except Exception as e:
            logger.error(f"Error saving results: {e}")
            raise
            
    def _prepare_results_for_serialization(self, results: Dict) -> Dict:
        """Prepare results dictionary for JSON serialization."""
        serializable = {}
        
        for key, value in results.items():
            if isinstance(value, (str, int, float, bool, list, dict)):
                serializable[key] = value
            elif isinstance(value, (np.integer, np.floating)):
                serializable[key] = value.item()
            elif isinstance(value, np.ndarray):
                serializable[key] = value.tolist()
            else:
                serializable[key] = str(value)
                
        return serializable
        
    def get_status(self) -> Dict:
        """Get current status of pipeline execution."""
        return {
            'timestamp': self.timestamp,
            'completed_steps': [k for k, v in self.results.items() if v],
            'output_directory': str(self.output_dir),
            'cache_directory': str(self.cache_dir)
        }

    def __del__(self):
        """Cleanup method to ensure HPC resources are released."""
        if self.hpc_manager:
            self.hpc_manager.stop_cluster()

def main():
    """
    Example usage of the SprintGenomePipeline.
    Demonstrates how to set up and run the complete analysis pipeline.
    """
    try:
        # Example paths - these would be replaced with actual paths
        vcf_path = "public/GFX0436892.filtered.snp.vcf.gz"
        output_dir = "output"
        config_path = "config/pipeline_config.json"
        
        logger.info("Initializing Sprint Genome Analysis Pipeline")
        
        # Initialize pipeline
        pipeline = SprintGenomePipeline(
            vcf_path=vcf_path,
            output_dir=output_dir,
            config_path=config_path,
            use_hpc=True  # Enable HPC for faster processing
        )
        
        # Run complete analysis
        results = pipeline.run_pipeline()
        
        # Get and display pipeline status
        status = pipeline.get_status()
        logger.info("Pipeline Status:")
        logger.info(f"Completed Steps: {status['completed_steps']}")
        logger.info(f"Results Directory: {status['output_directory']}")
        
        # Example of accessing specific results
        if results.get('genome_scoring'):
            logger.info(f"Total Sprint Score: {results['genome_scoring']['total_score']:.2f}")
            logger.info(f"Number of analyzed variants: {len(results['genome_scoring']['variant_scores'])}")
            
        if results.get('network_analysis'):
            logger.info(f"Network Size: {results['network_analysis']['summary']['num_nodes']} nodes, "
                       f"{results['network_analysis']['summary']['num_edges']} edges")
            
    except Exception as e:
        logger.error(f"Pipeline execution failed: {e}")
        raise

if __name__ == "__main__":
    main()
