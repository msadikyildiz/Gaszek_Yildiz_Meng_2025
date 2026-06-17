# Generate G matrix for epistasis analysis representing genotype mapping
import torch
from config import CACHE_DIR, get_mutation_options
from device_utils import get_device, clear_memory

def generate_G_matrix(mutation_options, device):
    """
    Generate matrix G that maps between different encodings of genotypes.
    Each mutation site contributes a block matrix based on number of variants.
    """
    G = torch.tensor([[1]], device=device, dtype=torch.float16)

    for option in mutation_options:
        if option == 2:
            B = torch.tensor([[1, 0],
                            [-1, 1]], device=device, dtype=torch.float16)
        elif option == 3:
            B = torch.tensor([[1, 0, 0],
                            [-1, 1, 0],
                            [-1, 0, 1]], device=device, dtype=torch.float16)
        else:  # option == 4
            B = torch.tensor([[1, 0, 0, 0],
                            [-1, 1, 0, 0],
                            [-1, 0, 1, 0],
                            [-1, 0, 0, 1]], device=device, dtype=torch.float16)
        
        G = torch.kron(B, G)  # Build matrix with Kronecker product
        del B
    
    return G

def main():
    """Generate and save the G matrix to disk"""
    device = get_device()
    mutation_opts = get_mutation_options()
    
    G = generate_G_matrix(mutation_opts, device)
    torch.save(G, CACHE_DIR / "G.pt")
    
    del G
    clear_memory(device)

if __name__ == "__main__":
    main()
