import os
 
def assert_cwd_is_bmeg():
    """Validates cwd is root dir of bmeg-etl."""
    assert 'transform' in [f for f in os.listdir('.') if os.path.isdir(f)]
    

    