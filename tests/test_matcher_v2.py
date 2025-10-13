"""
Tests for PubChem Matcher v2
=============================

Tests for ML-based matching with confidence scoring.
"""

import pytest
import json
from pathlib import Path

from matcher.matcher_v2 import MatcherV2, MatchResult
from matcher.confidence import MatchConfidenceEvaluator, Reliability
from matcher.similarity import MultiMetricSimilarity, SimilarityScore


# Test data
CLUSTER_GOOD = {
    'formula': 'C2H4O2',
    'atoms': ['C', 'C', 'O', 'O', 'H', 'H', 'H', 'H'],
    'bonds': [(0, 1), (0, 2), (1, 3), (0, 4), (0, 5), (1, 6), (1, 7)],
    'energy': -150.5,
    'positions': None
}

CLUSTER_FORMALDEHYDE = {
    'formula': 'CH2O',
    'atoms': ['C', 'H', 'H', 'O'],
    'bonds': [(0, 1), (0, 2), (0, 3, 2)],  # Double bond to O
    'energy': -114.2
}

CLUSTER_WATER = {
    'formula': 'H2O',
    'atoms': ['O', 'H', 'H'],
    'bonds': [(0, 1), (0, 2)],
    'energy': -76.4
}


class TestConfidenceEvaluator:
    """Test confidence evaluation"""
    
    def test_valence_check_valid(self):
        """Test valence checking with valid molecule"""
        evaluator = MatchConfidenceEvaluator()
        
        # Water: O-H-H (valence OK)
        cluster = {
            'atoms': ['O', 'H', 'H'],
            'bonds': [(0, 1, 1), (0, 2, 1)]
        }
        
        assert evaluator.check_valence(cluster) == True
    
    def test_valence_check_invalid(self):
        """Test valence checking with invalid molecule"""
        evaluator = MatchConfidenceEvaluator()
        
        # Carbon with 5 bonds (invalid!)
        cluster = {
            'atoms': ['C', 'H', 'H', 'H', 'H', 'H'],
            'bonds': [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1), (0, 5, 1)]
        }
        
        assert evaluator.check_valence(cluster) == False
    
    def test_charge_balance(self):
        """Test charge balance checking"""
        evaluator = MatchConfidenceEvaluator()
        
        # Neutral molecule
        cluster_neutral = {
            'atoms': [{'element': 'C', 'charge': 0}, {'element': 'O', 'charge': 0}],
            'bonds': []
        }
        assert evaluator.check_charge_balance(cluster_neutral) == True
        
        # High charge (should fail)
        cluster_charged = {
            'atoms': [
                {'element': 'C', 'charge': 2},
                {'element': 'O', 'charge': 2},
                {'element': 'N', 'charge': 2}
            ],
            'bonds': []
        }
        assert evaluator.check_charge_balance(cluster_charged) == False
    
    def test_bond_orders(self):
        """Test bond order checking"""
        evaluator = MatchConfidenceEvaluator()
        
        # Valid bond orders (1, 2, 3)
        cluster_valid = {
            'atoms': ['C', 'C', 'O'],
            'bonds': [(0, 1, 2), (1, 2, 1)]  # C=C-O
        }
        assert evaluator.check_bond_orders(cluster_valid) == True
        
        # Invalid bond order
        cluster_invalid = {
            'atoms': ['C', 'C'],
            'bonds': [(0, 1, 5)]  # Bond order 5!
        }
        assert evaluator.check_bond_orders(cluster_invalid) == False
    
    def test_evaluate_match_high_confidence(self):
        """Test match evaluation with high confidence"""
        evaluator = MatchConfidenceEvaluator()
        
        cluster = CLUSTER_GOOD
        match = {
            'cid': 757,
            'name': 'Glycolic acid',
            'molecular_formula': 'C2H4O2',
            'smiles': 'C(C(=O)O)O',
            'heavy_atom_count': 4,
            'similarity': 0.92
        }
        
        similarity = SimilarityScore(
            topology=0.95,
            fingerprint=0.90,
            energy=0.88,
            spectral=0.92,
            geometric=0.85
        )
        
        confidence = evaluator.evaluate_match(cluster, match, similarity)
        
        assert confidence.confidence_score >= 0.8
        assert confidence.reliability in [Reliability.HIGH, Reliability.MEDIUM]
        assert confidence.validation_status in ["PASS", "WARNING"]
        assert confidence.valence_check == True
        assert confidence.charge_balance == True
    
    def test_evaluate_match_low_confidence(self):
        """Test match evaluation with low confidence"""
        evaluator = MatchConfidenceEvaluator()
        
        cluster = CLUSTER_GOOD
        match = {
            'cid': 123,
            'name': 'Unknown',
            'molecular_formula': 'C3H6O',  # Different!
            'smiles': 'CCC=O',
            'heavy_atom_count': 5,
            'similarity': 0.35
        }
        
        similarity = SimilarityScore(
            topology=0.45,
            fingerprint=0.30,
            energy=0.40,
            spectral=0.25,
            geometric=0.35
        )
        
        confidence = evaluator.evaluate_match(cluster, match, similarity)
        
        assert confidence.confidence_score < 0.5
        assert confidence.reliability in [Reliability.LOW, Reliability.MEDIUM]
        assert len(confidence.warnings) > 0  # Should have warnings


class TestSimilarity:
    """Test multi-metric similarity"""
    
    def test_topology_similarity_identical(self):
        """Test topology similarity with identical molecules"""
        similarity = MultiMetricSimilarity()
        
        mol1 = {'formula': 'C2H4O2', 'bonds': [(0, 1), (1, 2)]}
        mol2 = {'formula': 'C2H4O2', 'bonds': [(0, 1), (1, 2)]}
        
        score = similarity.topology_similarity(mol1, mol2)
        assert score > 0.9  # Should be very similar
    
    def test_topology_similarity_different(self):
        """Test topology similarity with different molecules"""
        similarity = MultiMetricSimilarity()
        
        # Much more different molecules
        mol1 = {'formula': 'C2H4O2', 'bonds': [(0, 1), (1, 2), (2, 3), (3, 4)]}  # 4 bonds
        mol2 = {'formula': 'H2O', 'bonds': [(0, 1)]}  # 1 bond
        
        score = similarity.topology_similarity(mol1, mol2)
        assert score < 0.7  # Should be less similar
    
    def test_energy_similarity(self):
        """Test energy similarity"""
        similarity = MultiMetricSimilarity()
        
        # Similar energies
        mol1 = {'energy': -150.0, 'atoms': ['C', 'C', 'O', 'O']}
        mol2 = {'energy': -148.0, 'atoms': ['C', 'C', 'O', 'O']}
        
        score = similarity.energy_similarity(mol1, mol2)
        assert score > 0.8  # Should be similar
        
        # Different energies
        mol3 = {'energy': -100.0, 'atoms': ['C', 'C']}
        score2 = similarity.energy_similarity(mol1, mol3)
        assert score2 < score  # Should be less similar
    
    def test_compute_overall(self):
        """Test overall similarity computation"""
        similarity = MultiMetricSimilarity()
        
        mol1 = CLUSTER_GOOD
        mol2 = {
            'formula': 'C2H4O2',
            'atoms': ['C', 'C', 'O', 'O'],
            'bonds': [(0, 1), (0, 2), (1, 3)],
            'energy': -148.0,
            'smiles': 'C(C(=O)O)O'
        }
        
        score = similarity.compute(mol1, mol2, include_geometric=False)
        
        assert 0.0 <= score.overall <= 1.0
        assert 0.0 <= score.topology <= 1.0
        assert 0.0 <= score.fingerprint <= 1.0
        assert 0.0 <= score.energy <= 1.0


class TestMatcherV2:
    """Test MatcherV2 integration"""
    
    def test_init_without_ml(self):
        """Test initialization without ML classifier"""
        matcher = MatcherV2(
            classifier_model=None,
            use_ml_classifier=False
        )
        
        assert matcher.classifier is None
        assert matcher.use_ml == False
    
    def test_init_with_ml(self):
        """Test initialization with ML classifier"""
        # Skip if model doesn't exist
        model_path = 'data/atom_classifier.pkl'
        if not Path(model_path).exists():
            pytest.skip(f"Model not found: {model_path}")
        
        matcher = MatcherV2(
            classifier_model=model_path,
            use_ml_classifier=True
        )
        
        assert matcher.classifier is not None
        assert matcher.use_ml == True
    
    @pytest.mark.slow
    @pytest.mark.requires_network
    def test_match_cluster_water(self):
        """Test matching water molecule (requires network)"""
        matcher = MatcherV2(use_ml_classifier=False)  # Skip ML for speed
        
        result = matcher.match_cluster(CLUSTER_WATER, top_n=3)
        
        assert isinstance(result, MatchResult)
        assert result.cluster_formula == 'H2O'
        
        # Water should match successfully
        if result.success:
            assert result.pubchem_cid is not None
            assert 'water' in result.pubchem_name.lower() or result.pubchem_formula == 'H2O'
    
    @pytest.mark.slow
    @pytest.mark.requires_network
    def test_match_cluster_formaldehyde(self):
        """Test matching formaldehyde (requires network)"""
        matcher = MatcherV2(use_ml_classifier=False)
        
        result = matcher.match_cluster(CLUSTER_FORMALDEHYDE, top_n=3)
        
        assert isinstance(result, MatchResult)
        
        if result.success:
            assert result.pubchem_formula == 'CH2O'
            assert result.similarity_score.overall > 0.3
    
    def test_failed_result_creation(self):
        """Test creation of failed result"""
        matcher = MatcherV2(use_ml_classifier=False)
        
        result = matcher._failed_result(
            cluster=CLUSTER_GOOD,
            error_message="Test error"
        )
        
        assert result.success == False
        assert result.error_message == "Test error"
        assert result.pubchem_cid is None
        assert result.confidence.reliability == Reliability.INVALID
    
    def test_result_to_dict(self):
        """Test MatchResult serialization to dict"""
        matcher = MatcherV2(use_ml_classifier=False)
        
        result = matcher._failed_result(CLUSTER_GOOD, "Test")
        data = result.to_dict()
        
        assert isinstance(data, dict)
        assert 'cluster' in data
        assert 'pubchem' in data
        assert 'similarity' in data
        assert 'confidence' in data
        assert data['success'] == False
    
    def test_export_result(self, tmp_path):
        """Test exporting result to JSON"""
        matcher = MatcherV2(use_ml_classifier=False)
        
        result = matcher._failed_result(CLUSTER_GOOD, "Test")
        
        output_path = tmp_path / "test_result.json"
        matcher.export_result(result, str(output_path))
        
        assert output_path.exists()
        
        # Verify JSON is valid
        with open(output_path, 'r') as f:
            data = json.load(f)
        
        assert 'cluster' in data
        assert 'success' in data


# Run tests
if __name__ == "__main__":
    pytest.main([__file__, "-v", "-m", "not slow and not requires_network"])

