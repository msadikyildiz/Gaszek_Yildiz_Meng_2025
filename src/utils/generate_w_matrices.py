# Generate W matrices for different epistasis calculation methods
import torch
from config import CACHE_DIR
from device_utils import get_device, clear_memory

def calculate_W_linregress(V, X, S, H, device):
    """
    Calculate linear regression-based W matrix for epistasis calculation.
    Uses an optimized multi-step matrix multiplication approach with memory management.
    """
    temp = torch.matmul(S, H).to(torch.float16)
    del S, H
    clear_memory(device)
    
    temp2 = torch.matmul(X.T, temp).to(torch.float16)
    del temp, X
    clear_memory(device)
    
    result = torch.matmul(V, temp2).to(torch.float16)
    del temp2, V
    clear_memory(device)
    
    return result

def calculate_W_bioch(V, X, H):
    """
    Calculate biochemical W matrix for epistasis - more direct approach 
    requiring fewer matrices than linear regression method.
    """
    temp = torch.matmul(X.T, H).to(torch.float16)
    return torch.matmul(V, temp).to(torch.float16)

def main():
    """Generate and save W matrices for biochemical and linear regression epistasis models"""
    device = get_device()
    
    # Load necessary matrices and convert to float16
    print("Loading matrices...")
    V = torch.load(CACHE_DIR / "V.pt").to(torch.float16)
    X = torch.load(CACHE_DIR / "X.pt").to(torch.float16)
    H = torch.load(CACHE_DIR / "H.pt").to(torch.float16)
    
    # Calculate biochemical W first as it needs fewer matrices
    print("Calculating W_bioch...")
    W_bioch = calculate_W_bioch(V, X, H)
    torch.save(W_bioch, CACHE_DIR / "W_bioch.pt")
    del W_bioch
    clear_memory(device)
    
    # Calculate Q and S for linear regression
    print("Calculating Q and S...")
    Q = torch.eye(X.shape[0], device=device, dtype=torch.float16)
    S = torch.matmul(Q, Q.T).to(torch.float16)
    torch.save(S, CACHE_DIR / "S.pt")
    torch.save(Q, CACHE_DIR / "Q.pt")
    
    del Q
    clear_memory(device)

    # Calculate W_linregress
    print("Calculating W_linregress...")
    W_linregress = calculate_W_linregress(V, X, S, H, device)
    torch.save(W_linregress, CACHE_DIR / "W_linregress.pt")
    
    # Clean up
    del V, X, H, S, W_linregress
    clear_memory(device)

if __name__ == "__main__":
    main()
