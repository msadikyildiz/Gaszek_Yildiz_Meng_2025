import torch
from config import CACHE_DIR, get_mutation_options
from device_utils import get_device, clear_memory

def calculate_H_V(mutation_options, device):
    """
    Generate H and V matrices for epistasis analysis.
    H: Hadamard-like matrix representing genotype effects
    V: Matrix of coefficients for calculating epistatic effects
    """
    H = torch.tensor([[1.0]], dtype=torch.float16, device=device)
    V = torch.tensor([[1.0]], dtype=torch.float16, device=device)
    
    for option in mutation_options:
        if option == 2:
            B = torch.tensor([[1, 1], [1, -1]], device=device, dtype=torch.float16)
            C = torch.tensor([[1/2, 0], [0, -1]], device=device, dtype=torch.float16)
        elif option == 3:
            B = torch.tensor([[1, 1, 1], [1, -1, 0], [1, 0, -1]], device=device, dtype=torch.float16)
            C = torch.tensor([[1/3, 0, 0], [0, -1, 0], [0, 0, -1]], device=device, dtype=torch.float16)
        else:  # option == 4
            B = torch.tensor([[1, 1, 1, 1], [1, -1, 0, 0], [1, 0, -1, 0], [1, 0, 0, -1]], device=device, dtype=torch.float16)
            C = torch.tensor([[1/4, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]], device=device, dtype=torch.float16)
        
        H = torch.kron(H, B)  # Build H with Kronecker product
        V = torch.kron(V, C)  # Build V with Kronecker product
        del B, C
    
    return H, V

def calculate_W_ea(H, V):
    """Calculate ensemble-average epistasis matrix W_ea"""
    return torch.matmul(V, H)

def main():
    """Generate and save H, V matrices and ensemble-average W_ea matrix"""
    device = get_device()
    mutation_opts = get_mutation_options()
    
    H, V = calculate_H_V(mutation_opts, device)
    torch.save(H, CACHE_DIR / "H.pt")
    torch.save(V, CACHE_DIR / "V.pt")
    
    W_ea = calculate_W_ea(H, V)
    torch.save(W_ea, CACHE_DIR / "W_ea.pt")
    
    del H, V, W_ea
    clear_memory(device)

if __name__ == "__main__":
    main()