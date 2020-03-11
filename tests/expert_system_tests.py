
from common_util.expert_system import Expert_system
from nose.tools import *


class TestExpert_system(object):
    def test_loading(self):
        es = Expert_system()
        print(es.data)
        print(es.masses)
        m = es.masses
        m.sort()
        assert_equal(17.02654910101, m[0])
