
from mgf_checker.dataContainer import DataContainer


class MS_parser_info(DataContainer):
    def __init__(self, nQuery, dScore, stPeptide, stVarModsStr, stReadableVarMods, nMsmsid, stProtein_match, stFilename, nRank, stSeriesUsedStr):
        self.stReadableVarMods = stReadableVarMods
        self.nQuery = nQuery
        self.dScore = dScore
        self.stPeptide = stPeptide
        self.stVarModsStr = stVarModsStr
        self.nMsmsid = nMsmsid
        self.stProtein_match = stProtein_match
        self.stFilename = stFilename
        self.nRank = nRank
        self.stSeriesUsedStr = stSeriesUsedStr
