import requests
import pandas as pd
from typing import Dict, List, Optional, Union
from pathlib import Path
import logging
from dask import delayed, compute
import time
from ratelimit import limits, sleep_and_retry
from src.utils.hpc import HPCManager

logger = logging.getLogger(__name__)

class DatabaseIntegrator:
    """
    Integrates with various biological databases to enrich genetic analysis.
    Handles API rate limiting and caching of results.
    """
    
    def __init__(self, cache_dir: Optional[Path] = None):
        self.cache_dir = cache_dir or Path('cache')
        self.cache_dir.mkdir(exist_ok=True)
        self.hpc_manager = HPCManager()
        
        # API endpoints
        self.api_endpoints = {
            'string': 'https://string-db.org/api',
            'uniprot': 'https://rest.uniprot.org/uniprotkb',
            'reactome': 'https://reactome.org/ContentService/data',
            'kegg': 'https://rest.kegg.jp'
        }
        
    @sleep_and_retry
    @limits(calls=10, period=1)  # Rate limit: 10 calls per second
    def _make_api_request(self, url: str, params: Dict = None) -> Dict:
        """Make rate-limited API request with error handling."""
        try:
            response = requests.get(url, params=params)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            logger.error(f"API request failed: {e}")
            return {}
            
    def _get_cached_data(self, cache_key: str) -> Optional[pd.DataFrame]:
        """Retrieve cached data if available."""
        cache_file = self.cache_dir / f"{cache_key}.parquet"
        if cache_file.exists():
            try:
                return pd.read_parquet(cache_file)
            except Exception as e:
                logger.warning(f"Error reading cache: {e}")
        return None
        
    def _cache_data(self, data: pd.DataFrame, cache_key: str) -> None:
        """Cache data for future use."""
        try:
            cache_file = self.cache_dir / f"{cache_key}.parquet"
            data.to_parquet(cache_file)
        except Exception as e:
            logger.warning(f"Error caching data: {e}")
            
    @delayed
    def get_protein_interactions(self, gene: str) -> Dict:
        """Get protein-protein interactions from STRING database."""
        cache_key = f"string_{gene}"
        cached_data = self._get_cached_data(cache_key)
        
        if cached_data is not None:
            return cached_data.to_dict()
            
        url = f"{self.api_endpoints['string']}/network"
        params = {
            'identifier': gene,
            'species': 9606,  # Human
            'required_score': 700  # High confidence
        }
        
        data = self._make_api_request(url, params)
        if data:
            df = pd.DataFrame(data)
            self._cache_data(df, cache_key)
            return df.to_dict()
        return {}
        
    @delayed
    def get_pathway_data(self, gene: str) -> Dict:
        """Get pathway information from Reactome."""
        cache_key = f"reactome_{gene}"
        cached_data = self._get_cached_data(cache_key)
        
        if cached_data is not None:
            return cached_data.to_dict()
            
        url = f"{self.api_endpoints['reactome']}/pathways/entity/{gene}/contained"
        data = self._make_api_request(url)
        
        if data:
            df = pd.DataFrame(data)
            self._cache_data(df, cache_key)
            return df.to_dict()
        return {}
        
    def enrich_gene_data(self, genes: List[str]) -> Dict:
        """Enrich gene data using parallel API calls."""
        try:
            with self.hpc_manager.client:
                # Create tasks for parallel API calls
                ppi_tasks = []
                pathway_tasks = []
                
                for gene in genes:
                    ppi_task = self.hpc_manager.client.submit(
                        self.get_protein_interactions,
                        gene
                    )
                    pathway_task = self.hpc_manager.client.submit(
                        self.get_pathway_data,
                        gene
                    )
                    ppi_tasks.append(ppi_task)
                    pathway_tasks.append(pathway_task)
                
                # Gather results
                ppi_results = [task.result() for task in ppi_tasks]
                pathway_results = [task.result() for task in pathway_tasks]
                
                # Process results
                enriched_data = {
                    'protein_interactions': {
                        gene: ppi for gene, ppi in zip(genes, ppi_results)
                    },
                    'pathways': {
                        gene: pathway for gene, pathway in zip(genes, pathway_results)
                    }
                }
                
                return enriched_data
                
        except Exception as e:
            logger.error(f"Error in parallel data enrichment: {e}")
            raise
        
    def get_protein_metadata(self, uniprot_ids: List[str]) -> Dict:
        """Get detailed protein information from UniProt."""
        results = {}
        
        for uniprot_id in uniprot_ids:
            cache_key = f"uniprot_{uniprot_id}"
            cached_data = self._get_cached_data(cache_key)
            
            if cached_data is not None:
                results[uniprot_id] = cached_data.to_dict()
                continue
                
            url = f"{self.api_endpoints['uniprot']}/search"
            params = {
                'query': f'accession:{uniprot_id}',
                'format': 'json'
            }
            
            data = self._make_api_request(url, params)
            if data:
                df = pd.DataFrame(data['results'])
                self._cache_data(df, cache_key)
                results[uniprot_id] = df.to_dict()
                
        return results
        
    def get_metabolic_pathways(self, gene: str) -> Dict:
        """Get metabolic pathway information from KEGG."""
        cache_key = f"kegg_{gene}"
        cached_data = self._get_cached_data(cache_key)
        
        if cached_data is not None:
            return cached_data.to_dict()
            
        url = f"{self.api_endpoints['kegg']}/get/{gene}"
        data = self._make_api_request(url)
        
        if data:
            df = pd.DataFrame([data])
            self._cache_data(df, cache_key)
            return df.to_dict()
        return {}
