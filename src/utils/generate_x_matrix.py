import torch
from config import CACHE_DIR, get_mutation_options
from device_utils import get_device, clear_memory

def generate_X_matrix(mutation_options, device):
    """
    Generate matrix X for epistasis models, representing genotype encodings.
    Each row corresponds to a genotype, with columns representing mutation states.
    """
    X = torch.tensor([[1]], device=device, dtype=torch.float16)

    for option in mutation_options:
        if option == 2:
            D = torch.tensor([[1, 0],
                            [1, 1]], device=device, dtype=torch.float16)
        elif option == 3:
            D = torch.tensor([[1, 0, 0],
                            [1, 1, 0],
                            [1, 0, 1]], device=device, dtype=torch.float16)
        else:  # option == 4
            D = torch.tensor([[1, 0, 0, 0],
                            [1, 1, 0, 0],
                            [1, 0, 1, 0],
                            [1, 0, 0, 1]], device=device, dtype=torch.float16)
        
        X = torch.kron(X, D)  # Kronecker product to build combinatorial space
        del D
    
    return X

def main():
    """Generate and save the X matrix for all genotype combinations"""
    device = get_device()
    mutation_opts = get_mutation_options()
    
    X = generate_X_matrix(mutation_opts, device)
    torch.save(X, CACHE_DIR / "X.pt")
    
    del X
    clear_memory(device)

if __name__ == "__main__":
    main()