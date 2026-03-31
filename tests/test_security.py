import pytest
from main import app

@pytest.fixture
def client():
    app.config['TESTING'] = True
    with app.test_client() as client:
        yield client

def test_isoc_validation(client):
    """Test that /isoc route validates formula correctly."""
    # Valid formula
    response = client.post('/isoc', data={'expression': 'H2O', 'pts': '100', 'sigma': '0.1'})
    assert response.status_code == 200
    
    # Invalid formula (contains potentially malicious characters)
    response = client.post('/isoc', data={'expression': 'H2O<script>', 'pts': '100', 'sigma': '0.1'})
    assert response.status_code == 400

def test_m2f_validation(client):
    """Test that /m2f route validates mass correctly."""
    # Valid mass
    response = client.post('/m2f', data={'expression': '18.01', 'tol': '1'})
    assert response.status_code == 200
    
    # Invalid mass (contains non-numeric characters)
    response = client.post('/m2f', data={'expression': '18.01abc', 'tol': '1'})
    assert response.status_code == 400
    
    # Empty mass
    response = client.post('/m2f', data={'expression': '', 'tol': '1'})
    assert response.status_code == 400
