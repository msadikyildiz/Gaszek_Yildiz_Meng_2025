# Utility functions for GPU/CPU device management and memory cleanup
import torch
import gc

def get_device():
    """Detect and return the best available device (CUDA, MPS, or CPU)"""
    if torch.cuda.is_available():
        device = torch.device('cuda')
        torch.backends.cuda.matmul.allow_tf32 = True
        torch.backends.cudnn.allow_tf32 = True
        print("Using CUDA device.")
    elif getattr(torch.backends, "mps", None) and torch.backends.mps.is_available():
        device = torch.device('mps')
        print("Using Apple MPS device.")
    else:
        device = torch.device('cpu')
        print("Using CPU device.")
    return device

def clear_memory(device):
    """Free unused memory from Python garbage collector and device cache"""
    gc.collect()
    if device.type == 'cuda':
        torch.cuda.empty_cache()
        torch.cuda.synchronize()
    elif device.type == 'mps':
        torch.mps.empty_cache()