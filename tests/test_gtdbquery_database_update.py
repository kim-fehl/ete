import os
import pytest
import requests
from ete4 import GTDBTaxa, ETE_DATA_HOME
from ete4.gtdb_taxonomy import gtdbquery

DATABASE_PATH = ETE_DATA_HOME + '/gtdbtaxa.sqlite'
DEFAULT_GTDBTAXADUMP = ETE_DATA_HOME + '/gtdbdump.tar.gz'


@pytest.mark.network
def test_update_database():
    """Download and update the GTDB taxonomy database."""
    gtdb = GTDBTaxa()

    url = ('https://github.com/etetoolkit/ete-data/raw/main'
               '/gtdb_taxonomy/gtdb202/gtdb202dump.tar.gz')

    print(f'Downloading GTDB database release 202 to {DEFAULT_GTDBTAXADUMP} from {url}')

    with open(DEFAULT_GTDBTAXADUMP, 'wb') as f:
        f.write(requests.get(url).content)

    gtdb.update_taxonomy_database(DEFAULT_GTDBTAXADUMP)

    if not os.path.exists(DATABASE_PATH):
        gtdbquery.update_db(DATABASE_PATH)
