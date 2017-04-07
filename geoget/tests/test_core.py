from __future__ import division
import warnings
import unittest
import os
import shutil
import salem
from geoget.tests import is_download, is_slow, requires_credentials, cred
from geoget import core

# Setting for warnings
warnings.filterwarnings("once", category=DeprecationWarning)

# Globals
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_DIR = os.path.join(CURRENT_DIR, 'tmp_download')
if not os.path.exists(TEST_DIR):
    os.makedirs(TEST_DIR)


class TestFuncs(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass


class TestDataFiles(unittest.TestCase):

    def setUp(self):
        if os.path.exists(TEST_DIR):
            shutil.rmtree(TEST_DIR)
            os.makedirs(TEST_DIR)

    def tearDown(self):
        pass

    def test_download_sample_files(self):

        f = core.get_sample_file('OGGM/oggm-sample-data',
                                 'Hintereisferner.shp',
                                 TEST_DIR)
        self.assertTrue(os.path.exists(f))

        sh = salem.read_shapefile(f)
        self.assertTrue(hasattr(sh, 'geometry'))

    def test_srtmzone(self):

        z = core.srtm_zone(lon_ran=[-112, -112], lat_ran=[57, 57])
        self.assertTrue(len(z) == 1)
        self.assertEqual('14_01', z[0])

        z = core.srtm_zone(lon_ran=[-72, -73], lat_ran=[-52, -53])
        self.assertTrue(len(z) == 1)
        self.assertEqual('22_23', z[0])

        # Alps
        ref = sorted(['39_04', '38_03', '38_04', '39_03'])
        z = core.srtm_zone(lon_ran=[6, 14], lat_ran=[41, 48])
        self.assertTrue(len(z) == 4)
        self.assertEqual(ref, z)

    def test_asterzone(self):

        z, u = core.aster_zone(lon_ran=[137.5, 137.5], lat_ran=[-72.5, -72.5])
        self.assertTrue(len(z) == 1)
        self.assertTrue(len(u) == 1)
        self.assertEqual('S73E137', z[0])
        self.assertEqual('S75E135', u[0])

        z, u = core.aster_zone(lon_ran=[-95.5, -95.5],
                               lat_ran=[30.5, 30.5])
        self.assertTrue(len(z) == 1)
        self.assertTrue(len(u) == 1)
        self.assertEqual('N30W096', z[0])
        self.assertEqual('N30W100', u[0])

        z, u = core.aster_zone(lon_ran=[-96.5, -95.5],
                               lat_ran=[30.5, 30.5])
        self.assertTrue(len(z) == 2)
        self.assertTrue(len(u) == 2)
        self.assertEqual('N30W096', z[1])
        self.assertEqual('N30W100', u[1])
        self.assertEqual('N30W097', z[0])
        self.assertEqual('N30W100', u[0])

    def test_dem3_viewpano_zone(self):

        test_loc = {'ISL': [-25., -12., 63., 67.],  # Iceland
                    'SVALBARD': [10., 34., 76., 81.],
                    'JANMAYEN': [-10., -7., 70., 72.],
                    'FJ': [36., 66., 79., 82.],  # Franz Josef Land
                    'FAR': [-8., -6., 61., 63.],  # Faroer
                    'BEAR': [18., 20., 74., 75.],  # Bear Island
                    'SHL': [-3., 0., 60., 61.],  # Shetland
                    # Antarctica tiles as UTM zones, FILES ARE LARGE!!!!!
                    # '01-15': [-180., -91., -90, -60.],
                    # '16-30': [-91., -1., -90., -60.],
                    # '31-45': [-1., 89., -90., -60.],
                    # '46-60': [89., 189., -90., -60.],
                    # Greenland tiles
                    # 'GL-North': [-78., -11., 75., 84.],
                    # 'GL-West': [-68., -42., 64., 76.],
                    # 'GL-South': [-52., -40., 59., 64.],
                    # 'GL-East': [-42., -17., 64., 76.]
                    }
        # special names
        for key in test_loc:
            z = core.dem3_viewpano_zone([test_loc[key][0], test_loc[key][1]],
                                        [test_loc[key][2], test_loc[key][3]])
            self.assertTrue(len(z) == 1)

            self.assertEqual(key, z[0])

        # weird Antarctica tile
        # z = utils.dem3_viewpano_zone([-91., -90.], [-72., -68.])
        # self.assertTrue(len(z) == 1)
        # self.assertEqual('SR15', z[0])

        # normal tile
        z = core.dem3_viewpano_zone([-179., -178.], [65., 65.])
        self.assertTrue(len(z) == 1)
        self.assertEqual('Q01', z[0])

        # Alps
        ref = sorted(['K31', 'K32', 'K33', 'L31', 'L32',
                      'L33', 'M31', 'M32', 'M33'])
        z = core.dem3_viewpano_zone([6, 14], [41, 48])
        self.assertTrue(len(z) == 9)
        self.assertEqual(ref, z)

    @is_download
    def test_srtmdownload(self):

        # this zone does exist and file should be small enough for download
        zone = '68_11'
        fp = core.download_srtm_file(zone, TEST_DIR)
        self.assertTrue(os.path.exists(fp))
        fp = core.download_srtm_file(zone, TEST_DIR)
        self.assertTrue(os.path.exists(fp))

    @is_download
    def test_srtmdownloadfails(self):

        # this zone does not exist
        zone = '41_20'
        self.assertTrue(core.download_srtm_file(zone, TEST_DIR) is None)

    @is_download
    def test_iceland(self):
        fp, z = core.get_topo_file([-20, -20], [65, 65], TEST_DIR)
        self.assertTrue(os.path.exists(fp))

    @is_download
    def test_download_cru(self):

        of = core.get_cru_file(TEST_DIR, 'tmp')
        self.assertTrue(os.path.exists(of))

    @is_download
    @is_slow
    def test_download_rgi(self):

        of = core.get_rgi_data(TEST_DIR, version='5.0')
        of = os.path.join(of, '01_rgi50_Alaska', '01_rgi50_Alaska.shp')
        self.assertTrue(os.path.exists(of))

    @is_download
    def test_download_dem3_viewpano(self):

        # this zone does exist and file should be small enough for download
        zone = 'L32'
        fp = core.download_dem3_viewpano(zone, TEST_DIR)
        self.assertTrue(os.path.exists(fp))
        zone = 'U44'
        fp = core.download_dem3_viewpano(zone, TEST_DIR)
        self.assertTrue(os.path.exists(fp))

    @requires_credentials
    def test_get_postgresql_data(self):

        # take GLAMOS as test
        connect = cred['glamos']
        statement = "SELECT * FROM mass_balance.web_mass_balance_annual " \
                    "WHERE glacier_short_name = 'silvretta' AND xval = 1973;"
        df = core.get_postgresql_data(connect, statement)
        assert df.name.iloc[0] == 'Silvrettagletscher'

