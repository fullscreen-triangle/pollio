from typing import Callable, List, Any, Dict
from dask.distributed import Client, LocalCluster
import logging
from pathlib import Path
import os
import psutil
import numpy as np

logger = logging.getLogger(__name__)

class HPCManager:
    """
    Manages high-performance computing resources for genomic analysis.
    Provides utilities for parallel processing and memory management.
    """
    
    def __init__(self, n_workers: int = None, memory_limit: str = None):
        self.n_workers = n_workers or self._get_optimal_workers()
        self.memory_limit = memory_limit or self._get_memory_limit()
        self.client = None
        
    def _get_optimal_workers(self) -> int:
        """Determine optimal number of workers based on system resources."""
        cpu_count = os.cpu_count()
        return max(1, cpu_count - 1)  # Leave one CPU for system processes
        
    def _get_memory_limit(self) -> str:
        """Calculate memory limit per worker."""
        total_memory = psutil.virtual_memory().total
        worker_memory = int(total_memory * 0.8 / self.n_workers)  # 80% of available memory
        return f"{worker_memory}B"
        
    def start_cluster(self) -> None:
        """Start Dask local cluster with optimized settings."""
        try:
            cluster = LocalCluster(
                n_workers=self.n_workers,
                threads_per_worker=1,
                memory_limit=self.memory_limit
            )
            self.client = Client(cluster)
            logger.info(f"Started cluster with {self.n_workers} workers")
        except Exception as e:
            logger.error(f"Error starting cluster: {e}")
            raise
            
    def stop_cluster(self) -> None:
        """Shut down Dask cluster."""
        if self.client:
            self.client.close()
            self.client = None
            logger.info("Cluster shut down")
            
    def parallel_map(self, func: Callable, items: List[Any]) -> List[Any]:
        """Execute function in parallel across cluster."""
        if not self.client:
            self.start_cluster()
            
        try:
            futures = self.client.map(func, items)
            results = self.client.gather(futures)
            return results
        except Exception as e:
            logger.error(f"Error in parallel execution: {e}")
            raise
            
    def monitor_memory(self) -> Dict:
        """Monitor memory usage across workers."""
        if not self.client:
            return {}
            
        try:
            worker_info = self.client.scheduler_info()['workers']
            memory_usage = {
                worker: info['memory']
                for worker, info in worker_info.items()
            }
            return memory_usage
        except Exception as e:
            logger.error(f"Error monitoring memory: {e}")
            return {}
            
    @staticmethod
    def chunk_data(data: List[Any], chunk_size: int = None) -> List[List[Any]]:
        """Split data into chunks for parallel processing."""
        if chunk_size is None:
            chunk_size = max(1, len(data) // os.cpu_count())
            
        return [
            data[i:i + chunk_size]
            for i in range(0, len(data), chunk_size)
        ]
        
    def optimize_memory(self, data: Any) -> Any:
        """Optimize memory usage for large datasets."""
        if isinstance(data, (list, np.ndarray)):
            return self._optimize_array(data)
        return data
        
    def _optimize_array(self, arr: Any) -> Any:
        """Optimize array memory usage."""
        if isinstance(arr, list):
            arr = np.array(arr)
            
        if arr.dtype == np.float64:
            return arr.astype(np.float32)
        elif arr.dtype == np.int64:
            return arr.astype(np.int32)
        return arr 