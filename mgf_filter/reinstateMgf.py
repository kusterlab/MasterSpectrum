import sys
import blist


class mgfData(object):
    def __init__(self):
        pass

    def setTitle(self, title):
        self.title = str(title)

    def setPepmass(self, pepmass):
        self.pepmass = str(pepmass)

    def setCharge(self, charge):
        self.charge = str(charge)

    def setMZ(self, mz):
        self.mz = mz


class ReinstateMgf(object):
    def __init__(self, patched_mgf, original_mgf):
        self.patched_mgf = patched_mgf
        self.original_mgf = original_mgf
        self.patched_Data = blist.sorteddict()

    def writeMgf(self, path):
        with open(path, "wt") as file:
            file.write(self.firstLine)
            for key in self.patched_Data:
                pd = self.patched_Data[key]
                file.write("BEGIN IONS")
                file.write("\n")
                file.write("TITLE=" + pd.title)
                file.write("\n")
                file.write("PEPMASS=" + pd.pepmass)
                file.write("\n")
                file.write("CHARGE=" + pd.charge)
                file.write("\n")
                for d in pd.mz:
                    file.write(str(d) + '  ' + pd.mz[d][0] + '  ' + pd.mz[d][1])
                    file.write("\n")
                file.write("END IONS")
                file.write("\n\n")

    def compare2original(self):
        firstLine = True
        with open(self.original_mgf, "r") as file_to_read:
            for line in file_to_read:
                if firstLine:
                    firstLine = False
                    self.firstLine = line
                    continue

                if "BEGIN IONS" in line:
                    mD = mgfData()
                    mz_int_additional_data = blist.sorteddict()
                elif "END IONS" in line:
                    mD.setMZ(mz_int_additional_data)

                    for spectrum_keys in self.patched_Data[mD.title].mz.keys():
                        self.patched_Data[mD.title].mz[spectrum_keys] = mD.mz[spectrum_keys]
                else:
                    items = line.strip().split('=')
                    if items[0] == "TITLE":
                        mD.setTitle(items[1])
                    elif items[0] == "PEPMASS":
                        mD.setPepmass(items[1])
                    elif items[0] == "CHARGE":
                        mD.setCharge(items[1])
                    else:  # peaks
                        mz_int_additional = items[0].split()
                        if len(mz_int_additional) == 0:  # empty line
                            continue
                        else:
                            mz_int_additional_data[float(mz_int_additional[0])] = [mz_int_additional[1], mz_int_additional[2]]

    def loadPatched(self):
        with open(self.patched_mgf, "r") as file_to_read:
            for line in file_to_read:

                if "BEGIN IONS" in line:
                    mD = mgfData()
                    mz_int_additional_data = blist.sorteddict()
                elif "END IONS" in line:
                    mD.setMZ(mz_int_additional_data)
                    self.patched_Data[mD.title] = mD
                else:
                    items = line.strip().split('=')
                    if items[0] == "TITLE":
                        mD.setTitle(items[1])
                    elif items[0] == "PEPMASS":
                        mD.setPepmass(items[1])
                    elif items[0] == "CHARGE":
                        mD.setCharge(items[1])
                    else:  # peaks
                        mz_int_additional = items[0].split()
                        if len(mz_int_additional) == 0:  # empty line
                            continue
                        else:
                            mz_int_additional_data[float(mz_int_additional[0])] = [mz_int_additional[1]]


if __name__ == "__main__":
    rM = ReinstateMgf(sys.argv[1],
                      sys.argv[2])
    rM.loadPatched()
