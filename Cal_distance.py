import os
import re
from numpy import *


path = r"G:\msData\20200419\GST\DSS\output\reports"  # 搜索fasta 和 PDB 文件所在目录
XL_sites_list = ["K"]#, "D", 'E'] # 交联位点


###################Don't change the following line#####################
os.chdir(path)
# linksiteFile = "YLCao_DSS_BS3.csv" # site文件
for fl in os.listdir("./"):
    if fl == "pLink_summary.csv" or fl.endswith("v5.csv"):
        linksiteFile = fl


AA_dict = dict(
    HIS="H", MET="M", THR="T",
    PHE="F", PRO="P", SER="S",
    TRP="W", TYR="Y", VAL="V",
    M3L="m3K", GLY="G", ILE="I",
    ARG="R", LYS="K", LEU="L",
    ALA="A", CYS="C", ASN="N",
    GLN="Q", ASP="D", GLU='E')


def get_linked_site_inform(linked_site):
    pos_list = re.findall("\((\d*)\)", linked_site)
    position1 = pos_list[0]
    position2 = pos_list[1]
    p = linked_site.find("-")
    m = linked_site.find("(" + position1 + ")-")
    n = linked_site.find("(" + position2 + ")", p)
    protein1 = linked_site[:m].strip()
    protein2 = linked_site[p + 1:n].strip()
    return protein1, protein2, position1, position2


def find_maxvaule_poistion(matrix):
    max = matrix[0][0]
    for i in range(len(matrix[:, 0])):
        for j in range(len(matrix[0, :])):
            if matrix[i][j] > max:
                max = matrix[i][j]
                poistion = (i, j)
    return poistion, max


def longest_common_subsequence(lhs, rhs):
    l_len = len(lhs)
    r_len = len(rhs)
    matrix = zeros((l_len, r_len))
    for i in arange(l_len):
        for j in arange(r_len):
            if lhs[i] == rhs[j]:
                if i != 0 and j != 0:
                    matrix[i][j] = matrix[i - 1][j - 1] + 1
                else:
                    matrix[i][j] = 1
            elif i != 0 and j != 0:
                matrix[i][j] = max(matrix[i - 1][j], matrix[i][j - 1])
    return matrix


def longest_over_substrate(lhs, rhs):
    l_len = len(lhs)
    r_len = len(rhs)
    matrix = zeros((l_len, r_len))
    for i in arange(l_len):
        for j in arange(r_len):
            if lhs[i] == rhs[j]:
                if i != 0 and j != 0:
                    matrix[i][j] = matrix[i - 1][j - 1] + 1
                else:
                    matrix[i][j] = 1
    return matrix


def print_longest_over_substrate(lhs, rhs):
    matrix = longest_over_substrate(lhs, rhs)
    substrate_length = find_maxvaule_poistion(matrix)[1]
    return lhs[int(find_maxvaule_poistion(matrix)[0][0] - substrate_length) +
               1:int(find_maxvaule_poistion(matrix)[0][0]) + 1]


def print_longest_common_subsequence(lhs, rhs):
    matrix = longest_common_subsequence(lhs, rhs)
    l_len = len(lhs)
    r_len = len(rhs)
    i = l_len - 1
    j = r_len - 1
    rst = []
    while j > 0 and i > 0:
        if matrix[i][j] != matrix[i - 1][j] and matrix[i][j] != matrix[i][j
                                                                          - 1]:
            rst.insert(0, lhs[i])
            j -= 1
            i -= 1
        elif matrix[i][j] == matrix[i - 1][j]:
            i -= 1
        else:
            j -= 1
    if j != 0 and i == 0:
        if matrix[i][j - 1] != matrix[i][j]:
            rst.insert(0, lhs[i])
    elif j == 0 and 0:
        if matrix[i][j] != matrix[i - 1][j]:
            rst.insert(0, lhs[j])
    return rst


def find_xyz_coor(coordLine):
    x = float(coordLine[:8])
    y = float(coordLine[8:16])
    z = float(coordLine[16:24])
    return x,y,z


def cal_pdb_dis(chain_a, num1, chain_b, num2, strc_info_dic):
    coor1_line = strc_info_dic[chain_a, num1][-1]
    x1,y1,z1 = find_xyz_coor(coor1_line)
    coor2_line = strc_info_dic[chain_b, num2][-1]
    x2,y2,z2 = find_xyz_coor(coor2_line)
    
    sd = pow((x2 - x1),2) + pow((y2 - y1), 2) + pow((z2 - z1),2)
    return round(pow(sd, 0.5), 2)


def pretreatment_fasta(fasta):
    f = fasta
    fastaDic = {}
    i = 0
    while i < len(f):
        if f[i][0] != ">":
            i += 1
        else:
            nameLine = f[i]
            if " " in nameLine:
                namePro = nameLine.split(" ")[0][1:]
            else:
                namePro = nameLine[1:].strip()
            
            p = i + 1
            seq = ""
            while p < len(f):
                if f[p][0] == ">":
                    break
                else:
                    seq += f[p].strip()
                    p += 1
            fastaDic[namePro] = seq
            i = p
    
    return fastaDic


def trans2UpperSeq(line):
    seq_info = line[19:].rstrip()
    if len(seq_info) == 3:
        return AA_dict[seq_info]
    else:
        seq = ""
        seq_info_list = seq_info.split(" ")
        for aa in seq_info_list:
            seq += AA_dict[aa]
        return seq


def pretreatment_pdb(pdb):
    pdb_chain_To_seq_dic = {}
    for line in pdb:
        if line[:4] == "ATOM":
            break
        else:
            if line[:6] == "SEQRES":
                chain = line[11:13].strip()
                seq = trans2UpperSeq(line)
                if chain not in pdb_chain_To_seq_dic:
                    pdb_chain_To_seq_dic[chain] = seq
                else:
                    pdb_chain_To_seq_dic[chain] += seq
    
    pdb_chains = list(pdb_chain_To_seq_dic.keys())
    
    strc_info_dic = {}
    for line in pdb:
        if line[:4] == "ATOM":
            aa_num = int(line[22:26])
            chain =  line[21]
            atom_sym = line[13:15]
            aa = AA_dict[line[17:20]]
            if (chain, aa_num) not in strc_info_dic:  #元组需要加括弧
                strc_info_dic[chain, aa_num] = [aa]
                if atom_sym == "CA":
                    strc_info_dic[chain, aa_num].append(line[30:54])
            else:
                if atom_sym == "CA":
                    strc_info_dic[chain, aa_num].append(line[30:54])

    return strc_info_dic


# 根据chain序列逐一在fasta所有蛋白中匹配，找到匹配最好的，并返回匹配最好的蛋白以及序号差值
def judge_chain_from(chain_seq, fastaDic, real_index_list):
    score_dic = {}
    for prot in fastaDic:
        overlap = print_longest_over_substrate(fastaDic[prot], chain_seq)
        score = len(overlap)/len(chain_seq)
        score_dic[prot] = [score, overlap]
    new_list = sorted(list(score_dic.items()), key = lambda x:x[1][0])
    prot_max, (max_score, max_overlap) = new_list[-1]

    if max_score < 0.7:
        print("warining the overlap is {:.2f}".format(max_score))
    
    index_in_fasta = fastaDic[prot_max].index(max_overlap)
    index_in_pdb_seq = chain_seq.index(max_overlap)
    real_index_in_Pdb = real_index_list[index_in_pdb_seq]
    delta = index_in_fasta + 1 - real_index_in_Pdb
    return prot_max, delta


def get_fasta_pdb_info(fasta, pdb):
    fastaDic = pretreatment_fasta(fasta)
    strc_info_dic = pretreatment_pdb(pdb)
    key_list = sorted(strc_info_dic.keys(), key = lambda x:x[1])
    structuredChain_To_seq = {}
    chain_numList4seqInpdb = {}
    for (chain, aa_num) in key_list:
        if chain not in structuredChain_To_seq:
            chain_numList4seqInpdb[chain] = [aa_num]
            structuredChain_To_seq[chain] = strc_info_dic[chain, aa_num][0]
        else:
            chain_numList4seqInpdb[chain].append(aa_num)
            structuredChain_To_seq[chain] += strc_info_dic[chain, aa_num][0]
    
    chain2protName_dic = {} # chain与蛋白名称的对应关系
    pos_delta_pdb2fasta = {} #chain里fasta序号与pdb序号中的差值
    for chain in structuredChain_To_seq:
        real_index_list = chain_numList4seqInpdb[chain]
        chain_from, delta = judge_chain_from(structuredChain_To_seq[chain], fastaDic, real_index_list)
        chain2protName_dic[chain] = chain_from
        pos_delta_pdb2fasta[chain] = delta

    return chain2protName_dic, pos_delta_pdb2fasta


def find_chain_pos(protein, position, chain2protName_dic, pos_delta_pdb2fasta, strc_info_dic):
    prot = protein
    pos = position
    site_list = []
    for chain in chain2protName_dic:
        if chain2protName_dic[chain] == prot:
            Cor_Pos1 = int(pos) - pos_delta_pdb2fasta[chain]
            if (chain, Cor_Pos1) in strc_info_dic:
                match = strc_info_dic[chain, Cor_Pos1]
                if len(match) == 2:
                    aa, coord = match
                    if aa in XL_sites_list:
                        site_list.append((chain, Cor_Pos1))
    return site_list


# 对于给出的cross-link pair给出其最短距离
def get_pdb_distance(cross_link_pair, chain2protName_dic,
                     pos_delta_pdb2fasta, strc_info_dic, pdb_name):
    ptn1, ptn2, pos1, pos2 = get_linked_site_inform(cross_link_pair)
    site_1_chain = find_chain_pos(ptn1, pos1, chain2protName_dic, pos_delta_pdb2fasta, strc_info_dic)
    site_2_chain = find_chain_pos(ptn2, pos2, chain2protName_dic, pos_delta_pdb2fasta, strc_info_dic)
    
    print(site_1_chain, site_2_chain)
    wlist = [cross_link_pair]
    if len(site_1_chain) and len(site_2_chain):
        all_distance = []
        all_dist_dic = {}
        for i in range(len(site_1_chain)):
            chain1, site1 = site_1_chain[i]
            for j in range(len(site_2_chain)):            
                chain2, site2 = site_2_chain[j]

                dis_pdb = cal_pdb_dis(chain1, site1, chain2, site2, strc_info_dic)
                obj_name = "dist " + chain1 + str(pos1) + "_" + chain2 + str(pos2)
                selec1 = " /" + pdb_name + "//" + chain1 + "/" + str(site1) + "/"+"CA"
                selec2 = " /" + pdb_name + "//" + chain2 + "/" + str(site2) + "/"+"CA"
                pym_cmd = obj_name + "," + selec1 + "," + selec2
                all_dist_dic[site_1_chain[i], site_2_chain[j]] = [str(dis_pdb), pym_cmd]
                all_distance.append(dis_pdb)
                wlist.extend([dis_pdb, pym_cmd])
        
        if len(all_distance) == 1:
            min_dis = round(min(all_distance), 2)
        else:
            while min(all_distance) == 0.0:
                all_distance.remove(0.0)
            min_dis = round(min(all_distance), 2)
        wlist.insert(1, min_dis)
    else:
        wlist.append("no stru")
    
    return wlist


def main():
    B = open("xlink_distance.txt", 'w')
    B.write("link-site\ttotal_spec\tbest_Evaule\tbest_svm\tmin_distance\n")
    file_list = os.listdir(os.getcwd())
    for fl in file_list:
        if fl[-6:] == ".fasta":
            fasta = open(fl, 'r').readlines()
            print("The fasta file is " + fl)
        elif fl[-4:] == ".pdb":
            pdb_name = fl[:-4]
            pdb = open(fl, 'r').readlines()
            print("The pdb file is " + pdb_name)
        else:
            continue
    
    chain2protName_dic, pos_delta_pdb2fasta = get_fasta_pdb_info(fasta, pdb)
    strc_info_dic = pretreatment_pdb(pdb)
    print(chain2protName_dic)
    print(pos_delta_pdb2fasta)

    f = open(linksiteFile, 'r').readlines()
    
    for line in f[1:]:   #important
        line_list = line.rstrip("\n").split(",")
        # print(line_list)
        if line_list[0].isdigit():
            break
        else:
            link_pair = line_list[0]
            if "(" in link_pair and ")" in link_pair:
                wlist = line_list[0:4]
                wlist.extend(get_pdb_distance(link_pair, chain2protName_dic, pos_delta_pdb2fasta, strc_info_dic, pdb_name)[1:])
                # wlist.extend()
                B.write("\t".join([str(ele) if type(ele) != str else ele for ele in wlist])+"\n")
        
    B.close()


if __name__ == "__main__":
    main()
    print("Done")
