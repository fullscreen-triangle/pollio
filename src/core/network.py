import networkx as nx
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Set, Tuple
from pathlib import Path
import logging
from dask import delayed, compute
from src.utils.hpc import HPCManager

logger = logging.getLogger(__name__)

class NetworkBuilder:
    """
    Builds and analyzes biological networks related to sprint performance genes.
    Integrates protein-protein interactions, metabolic pathways, and regulatory relationships.
    """
    
    def __init__(self, genome_data: Dict, config_path: Optional[Path] = None):
        self.genome_data = genome_data
        self.network = nx.Graph()
        self.config = self._load_config(config_path) if config_path else {}
        self.ppi_network = nx.Graph()
        self.metabolic_network = nx.Graph()
        self.regulatory_network = nx.DiGraph()
        self.hpc_manager = HPCManager()
        
    def _load_config(self, config_path: Path) -> Dict:
        """Load network configuration parameters."""
        try:
            return pd.read_json(config_path).to_dict()
        except Exception as e:
            logger.error(f"Error loading network config: {e}")
            return {}
            
    @delayed
    def _build_ppi_subnetwork(self, gene: str) -> Set[Tuple[str, str]]:
        """Build protein-protein interaction subnetwork for a gene."""
        try:
            # This would typically make API calls to STRING or BioGRID
            # For now, returning placeholder edges
            edges = set()
            if gene in self.genome_data['variants']:
                # Add first-degree interactions
                edges.add((gene, f"{gene}_interactor1"))
                edges.add((gene, f"{gene}_interactor2"))
            return edges
        except Exception as e:
            logger.error(f"Error building PPI network for {gene}: {e}")
            return set()
            
    def build_network(self) -> nx.Graph:
        """Build network using parallel processing for subnetworks."""
        try:
            with self.hpc_manager.client:
                # Create tasks for parallel network construction
                network_tasks = []
                for gene in self.genome_data['variants'].keys():
                    task = self.hpc_manager.client.submit(
                        self._build_ppi_subnetwork,
                        gene
                    )
                    network_tasks.append(task)
                
                # Gather results
                subnetworks = [task.result() for task in network_tasks]
                
                # Merge subnetworks
                self.network = nx.compose_all([
                    subnetwork for subnetwork in subnetworks if subnetwork
                ])
                
                return self.network
                
        except Exception as e:
            logger.error(f"Error in parallel network construction: {e}")
            raise
            
    def _add_node_attributes(self) -> None:
        """Add relevant attributes to network nodes."""
        for node in self.network.nodes():
            if node in self.genome_data['variants']:
                variant_data = self.genome_data['variants'][node]
                nx.set_node_attributes(
                    self.network,
                    {node: {
                        'type': 'gene',
                        'score': variant_data.get('score', 0.0),
                        'beneficial': variant_data.get('beneficial', False)
                    }}
                )
                
    def analyze_network(self) -> Dict:
        """Perform network analysis to identify key features."""
        if not self.network:
            raise ValueError("Network not built. Call build_network first.")
            
        analysis_results = {
            'centrality': self._calculate_centrality(),
            'communities': self._detect_communities(),
            'pathways': self._analyze_pathways(),
            'summary': self._generate_network_summary()
        }
        
        return analysis_results
        
    def _calculate_centrality(self) -> Dict:
        """Calculate various centrality metrics for network nodes."""
        return {
            'degree': nx.degree_centrality(self.network),
            'betweenness': nx.betweenness_centrality(self.network),
            'eigenvector': nx.eigenvector_centrality_numpy(self.network)
        }
        
    def _detect_communities(self) -> List[Set]:
        """Detect communities in the network."""
        return list(nx.community.greedy_modularity_communities(self.network))
        
    def _analyze_pathways(self) -> List[Dict]:
        """Analyze significant pathways in the network."""
        pathways = []
        beneficial_genes = {
            node for node in self.network.nodes()
            if self.network.nodes[node].get('beneficial', False)
        }
        
        # Find paths between beneficial genes
        for source in beneficial_genes:
            for target in beneficial_genes:
                if source != target:
                    try:
                        path = nx.shortest_path(self.network, source, target)
                        pathways.append({
                            'source': source,
                            'target': target,
                            'path': path,
                            'length': len(path)
                        })
                    except nx.NetworkXNoPath:
                        continue
                        
        return pathways
        
    def _generate_network_summary(self) -> Dict:
        """Generate summary statistics for the network."""
        return {
            'num_nodes': self.network.number_of_nodes(),
            'num_edges': self.network.number_of_edges(),
            'density': nx.density(self.network),
            'avg_clustering': nx.average_clustering(self.network),
            'avg_shortest_path': nx.average_shortest_path_length(self.network)
        }
        
    def visualize_network(self, output_path: Path) -> None:
        """Generate network visualization."""
        try:
            import matplotlib.pyplot as plt
            
            plt.figure(figsize=(12, 12))
            pos = nx.spring_layout(self.network)
            
            # Draw nodes
            nx.draw_networkx_nodes(
                self.network, pos,
                node_color=[
                    'g' if self.network.nodes[node].get('beneficial', False) else 'r'
                    for node in self.network.nodes()
                ],
                node_size=500
            )
            
            # Draw edges
            nx.draw_networkx_edges(self.network, pos)
            
            # Add labels
            nx.draw_networkx_labels(self.network, pos)
            
            plt.title("Sprint Performance Genetic Network")
            plt.axis('off')
            
            plt.savefig(output_path / 'network_visualization.png')
            plt.close()
            
        except Exception as e:
            logger.error(f"Error visualizing network: {e}")
