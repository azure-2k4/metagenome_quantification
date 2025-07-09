#!/usr/bin/env python3
"""
Unit tests for quantification methods.
"""
import numpy as np
import pytest
from pathlib import Path
import sys

# Add scripts directory to Python path
sys.path.append(str(Path(__file__).parent.parent / 'scripts'))

from quantify_rpkm import calculate_rpkm, normalize_abundances as rpkm_normalize
from quantify_poisson import poisson_loglik
from quantify_nnls import build_presence_matrix, normalize_abundances as nnls_normalize

def test_rpkm_calculation():
    """Test RPKM calculation."""
    # Test case: 100 reads mapped to 1000bp sequence with 1M total reads
    rpkm = calculate_rpkm(100, 1000, 1000000)
    expected = 100
    assert np.isclose(rpkm, expected)
    
    # Test edge cases
    assert calculate_rpkm(0, 1000, 1000000) == 0
    with pytest.raises(ZeroDivisionError):
        calculate_rpkm(100, 0, 1000000)

def test_rpkm_normalization():
    """Test abundance normalization."""
    # Test case: Convert RPKM values to proportions
    rpkm_values = {'A': 100, 'B': 200, 'C': 700}
    normalized = rpkm_normalize(rpkm_values)
    
    assert np.isclose(sum(normalized.values()), 1.0)
    assert normalized['C'] > normalized['B'] > normalized['A']

def test_poisson_likelihood():
    """Test Poisson log-likelihood calculation."""
    # Simple test case with 2 genomes and 3 unitigs
    abundances = np.array([0.3, 0.7])
    read_counts = np.array([10, 20, 30])
    presence_matrix = np.array([
        [1, 0],
        [0, 1],
        [1, 1]
    ])
    unitig_lengths = np.array([100, 100, 100])
    total_reads = 60
    
    loglik = poisson_loglik(
        abundances,
        read_counts,
        presence_matrix,
        unitig_lengths,
        total_reads
    )
    
    # Log-likelihood should be finite and negative
    assert np.isfinite(loglik)
    assert loglik < 0

def test_nnls_matrix_construction():
    """Test construction of presence matrix for NNLS."""
    unitigs = ['u1', 'u2', 'u3']
    colors = {
        'u1': [1],
        'u2': [2],
        'u3': [1, 2]
    }
    
    A, genomes = build_presence_matrix(unitigs, colors)
    
    assert A.shape == (3, 2)  # 3 unitigs, 2 genomes
    assert np.array_equal(A[0], [1, 0])  # u1 belongs to genome 1
    assert np.array_equal(A[1], [0, 1])  # u2 belongs to genome 2
    assert np.array_equal(A[2], [1, 1])  # u3 belongs to both

def test_nnls_normalization():
    """Test normalization of NNLS abundances."""
    # Test case
    abundances = np.array([2, 4, 6])
    normalized = nnls_normalize(abundances)
    
    assert np.isclose(sum(normalized), 1.0)
    assert np.allclose(normalized, [0.166667, 0.333333, 0.5])

if __name__ == '__main__':
    pytest.main([__file__])
