import os
import re

wk_dir = r"G:\msData\20200419\GST\DSS\output\reports"
os.chdir(wk_dir)

def get_input_file():
    for fl in os.listdir(wk_dir):
        if fl.endswith("v5.csv") or fl == "pLink_summary.csv":
            return fl
    print("Wrong")

input_file = get_input_file()

def get_linked_site_inform(linked_site):
    pos_list = re.findall("\((\d*)\)", linked_site)
    position1 = pos_list[0]
    position2 = pos_list[1]
    inf1, inf2 = linked_site.split(")-")
    protein1 = inf1[:inf1.find("(" + position1)]
    protein2 = inf2[:inf2.find("("+position2)]
    # p = linked_site.find(")-")
    # m = linked_site.find("(" + position1 + ")-")
    # n = linked_site.find("(" + position2 + ")", p)
    # protein1 = linked_site[:m].strip()
    # protein2 = linked_site[p + 2:n].strip()
    if protein1 != protein2:
        link_type = "Inter"
    else:
        if abs(int(position1)-int(position2)) < 5:
            link_type = "Inter"
        else:
            link_type = "Intra" 
    return protein1, protein2, position1, position2, link_type


def main():
    file_list = os.listdir(os.getcwd())
    for file in file_list:
        if input_file == file:
            f = open(file).readlines()
            rep_name = file[:-4]+"_xinet.csv"
            b = open(rep_name, "w")
            b.write(",".join(["Score", "Protein1", "Protein2", "LinkPos1", "LinkPos2"]))
            b.write("\n")

            for line in f[1:]:
                line_list = line.rstrip("\n").split(",")
                site = line_list[0]
                spectra = line_list[1]
                write_list = [spectra]
                if site.isdigit():
                    break
                else:                
                    if "/" in site:
                        continue
                    else:
                        if "Molecular" in site:
                            continue
                        else:
                            site_extrac = get_linked_site_inform(site)
                            if site_extrac[-1] == "Inter" or site_extrac[-1] == "Intra" :
                                write_list.extend(site_extrac[:-1])
                                b.write(",".join(write_list))
                                print(write_list)
                                b.write("\n")
            b.close()
            print("Done")
        else:
            continue


if __name__ == "__main__":
    main()
