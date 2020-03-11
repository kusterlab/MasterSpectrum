
import blist
import os


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


class Apl2Mgf(object):
    def __init__(self, apl_folder, output_mgf):
        self.apl_folder = apl_folder
        self.output_mgf = output_mgf
        self.apl_data = blist.sorteddict()

    def writeMgf(self):
        with open(self.output_mgf, "wt") as file:
            for key in self.apl_data:
                pd = self.apl_data[key]
                file.write("BEGIN IONS")
                file.write("\n")
                file.write("TITLE=" + pd.title)
                file.write("\n")
                file.write("PEPMASS=" + pd.pepmass)
                file.write("\n")
                file.write("CHARGE=" + pd.charge + "+")
                file.write("\n")
                for d in pd.mz:
                    file.write(str(d) + '  ' + str(float(pd.mz[d][0])))
                    file.write("\n")
                file.write("END IONS")
                file.write("\n\n")

    def load_apls(self):
        for f in self.apl_files:
            file = self.apl_folder + f
            print(file)
            with open(file, "r") as file_to_read:
                for line in file_to_read:
                    line = line.strip()
                    if "peaklist start" in line:
                        mD = mgfData()
                        mz_int_additional_data = blist.sorteddict()
                    elif "peaklist end" in line:
                        mD.setMZ(mz_int_additional_data)
                        self.apl_data[mD.title + "__tobi__" + mD.charge] = mD
                    else:
                        items = line.strip().split('=')
                        if items[0] == 'mz':
                            mD.setPepmass(items[1])
                        elif items[0] == "header":
                            mD.setTitle(items[1])
                        elif items[0] == "charge":
                            mD.setCharge(items[1])
                        elif items[0] == "fragmentation":
                            continue
                        else:  # peaks
                            mz_int_additional = items[0].split("\t")
                            if len(mz_int_additional) == 1:  # empty line
                                continue
                            else:
                                mz_int_additional_data[float(mz_int_additional[0])] = [mz_int_additional[1]]
            print(len(self.apl_data))

    def find_all_apl(self):
        self.apl_files = [f for f in os.listdir(self.apl_folder) if f.endswith('.apl')]
        print(self.apl_files)
