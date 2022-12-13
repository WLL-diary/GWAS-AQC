import sys, os, subprocess, statistics
import requests
from flask import Flask, jsonify, request
import time
import requests
from flask import Flask, jsonify, request
time_start=time.time()
application = Flask(__name__)
fp = open("C:/Users/WLL/PycharmProjects/pythonProject/final_QC.log")
for num, line in enumerate(fp):
        if num == 6:
            a=line
            n0=a[4]+a[5]+a[6]
        if num == 7:
            a = line
            n1 = a[5] + a[6] + a[7]
        if num == 8:
            a = line
            n2 = a[5] + a[6] + a[7]
        if num == 9:
            a = line
            n3 = a[4] + a[5] + a[6]+a[7]+a[8]+a[9]
        if num==30:
            a=line
            n5=a[11]+a[12]+a[13]+a[14]+a[15]+a[16]+a[17]
        if num==40:
            a=line
            n4=a[5]+a[6]+a[7]

MAF=n0
MIND=n1
GENO=n2
HWE=n3
Individual=n4
variants=n5
#上传至区块链的内容
 # HTTP Transfer over Blockchain
Post_Transaction = 'http://localhost:5000/transactions/new'
Mine_Transaction = 'http://localhost:5000/mine'
#
#
def post_blockchain_data(MAF,MIND,GENO,HWE,Individual,variants):
    """
    Creates a JSON file format payload to send over the HTTP network and complete addition of a valid block to chain
    创建一个 JSON 文件格式的有效负载以通过 HTTP 网络发送并完成将有效块添加到链中
    """
    n0_text = str(MAF)
    n1_text = str(MIND)
    n2_text = str(GENO)
    n3_text = str(HWE)
    n4_text=str(Individual)
    n5_text=str(variants)
    payload ="{\n \"MAF\":\""+n0_text+"\",\"MIND\":\""+n1_text+"\",\"GENO\":\""+n2_text+"\",\"HWE\":\""+n3_text+"\",\"Individual\":\""+n4_text+"\",\"variants\":\""+n5_text+"\"}"
    print(payload)
    headers = {
        'Content-Type': "application/json",
    }
    response = requests.request("POST", Post_Transaction, data= payload, headers=headers)
    print(response.text)
    return

pfile=open(sys.argv[1], 'r')

params = {}
for line in pfile:
  paramname, paramvalue = line.split()
  params[paramname[1:]]=paramvalue

bfile=params['WORK']+'/'+params['INPUT'].split('/')[-1]
logfile=open(params['WORK']+'/final_QC.log', 'w')
logfile.write('Input parameters: '+sys.argv[1]+'\n')
for p in params:
  logfile.write(p+'\t'+params[p]+'\n')
logfile.write('************************************\n\n')

#CHECK SEX
#性别检查

if params['SEX'] == 'yes':
  os.system(params['PLINKPATH']+' --bfile '+params['INPUT']+' --check-sex --pheno '+params['PFILE']+' --pheno-name '+params['PHENO']+' --out '+bfile+'_sexcheck')
  #number of samples
  proc = subprocess.Popen(['findstr', "people", bfile+'_sexcheck.log'], stdout=subprocess.PIPE)
  tmp, err = proc.communicate()
  n = int(tmp.split()[0])
  logfile.write('n = '+str(n)+'\n')

  #number of variants
  #变体数量
  proc = subprocess.Popen(['findstr', "variants", bfile+'_sexcheck.log'], stdout=subprocess.PIPE)
  tmp, err = proc.communicate()
  v = int(tmp.split()[0])
  logfile.write('variants = '+str(v)+'\n\n')

  #problems with sex check
  #有问题的性别
  proc = subprocess.Popen(['findstr', "PROBLEM", bfile+'_sexcheck.sexcheck'], stdout=subprocess.PIPE)
  tmp, err = proc.communicate()

  tmp = tmp.decode('ascii')
  tmp_split=tmp.split()
  print(tmp)
  print(tmp_split)
  nproblems=int(len(tmp_split)/6)

  #remove problems
  #删除有问题的个体
  logfile.write('Check Sex Problems: '+str(nproblems)+'\n')
  if nproblems > 0:
    rfile=open(bfile+'_remove.txt', 'w')
    for k in range(nproblems):
      rfile.write(tmp_split[k*6]+'\t'+tmp_split[1+k*6]+'\n')
    rfile.close()
    os.system(params['PLINKPATH']+' --bfile '+params['INPUT']+' --remove '+bfile+'_remove.txt --make-bed --out '+bfile+'_sex')
    params['INPUT'] = bfile+'_sex'
    bfile=bfile+'_sex'
    proc = subprocess.Popen(['findstr', "remove:", bfile+'.log'], stdout=subprocess.PIPE)
    tmp = proc.stdout.read()
    n = int(tmp.split()[1])
    logfile.write('n =  '+str(n)+'\n\n')

#investigate missingness
#缺失质控
os.system(params['PLINKPATH']+' --bfile '+params['INPUT']+' --missing --pheno '+params['PFILE']+' --pheno-name '+params['PHENO']+'  --out '+bfile+'_miss')
os.system(params['RPATH']+' '+params['SCRIPTS']+'/missing.R '+bfile+'_miss.imiss '+bfile+'_miss.lmiss '+params['WORK']+ ' ')
#os.system('mv Rplots.pdf '+bfile+'_missing.pdf')
imiss = bfile+'_miss.imiss'

#if n and v not obtained during sex check, do here
if params['SEX'] == 'no':
  #number of samples
  #样本数量
  proc = subprocess.Popen(['findstr', "people", bfile+'_miss.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  n = int(tmp.split()[0])
  logfile.write('n = '+str(n)+'\n')

  #number of variants
  #变体数量
  proc = subprocess.Popen(['findstr', "variants", bfile+'_miss.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  v = int(tmp.split()[0])
  logfile.write('variants = '+str(v)+'\n\n')

#apply thresholds
#应用阈值
if params['MAF'] == 'na':
  params['MAF'] = str(10.0/n)
  logfile.write('MAF = 10/n = '+params['MAF']+'\n')

if 'MIND1' not in params:
  os.system(params['PLINKPATH']+' --bfile '+params['INPUT']+' --geno '+params['GENO']+' --maf '+params['MAF']+' --mind '+params['MIND']+' --make-bed --out '+bfile+'_thresh')
  #now bfile and input should always be same so just use bfile
  bfile=bfile+'_thresh'
  #report on QC
  proc = subprocess.Popen(['findstr', "mind)", bfile+'.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  mind = int(tmp.split()[0])
  logfile.write('Samples removed after --mind = '+str(mind)+'\n')

  if (float(params['GENO']) < 1):
    proc = subprocess.Popen(['findstr', "geno)", bfile+'.log'], stdout=subprocess.PIPE)
    tmp = proc.stdout.read()
    geno = int(tmp.split()[0])
  else:
    geno = 0
  logfile.write('Variants removed after --geno = '+str(geno)+'\n')
  proc = subprocess.Popen(['findstr "maf)"', bfile+'.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  maf = int(tmp.split()[0])
  logfile.write('Variants removed after --maf = '+str(maf)+'\n')
  proc = subprocess.Popen(['findstr', 'pass', bfile+'.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  oldv = v
  oldn = n
  v = int(tmp.split()[0])
  n = int(tmp.split()[3])
  logfile.write('\nn = '+str(n)+'\nvariants = '+str(v)+'\n')

  #check
  if oldv != (v+maf+geno):
    logfile.write('Check variant QC\n')
  if oldn != (n+mind):
    logfile.write('Check sample QC\n')

#2-tier MIND theshold
else:
#first MIND threshold
#使用--mind删除个体
  os.system(params['PLINKPATH']+' --bfile '+params['INPUT']+' --mind '+params['MIND1']+' --make-bed --out '+bfile+'_thresh')
  #now bfile and input should always be same sp just use bfile
  bfile=bfile+'_thresh'
  #report on QC
  proc = subprocess.Popen(['findstr', 'mind)', bfile+'.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  mind = int(tmp.split()[0])
  logfile.write('Samples removed after --mind(1) = '+str(mind)+'\n')
  proc = subprocess.Popen(['findstr', 'pass', bfile+'.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  oldn = n
  n = int(tmp.split()[3])
  logfile.write('\nn = '+str(n)+'\n')
  #check
  if oldn != (n+mind):
    logfile.write('Check sample QC\n')
#other thresholds
#使用--geno删除SNP
#使用--maf删除次等位基因
  os.system(params['PLINKPATH']+' --bfile '+bfile+' --geno '+params['GENO']+' --maf '+params['MAF']+' --make-bed --out '+bfile+'_thresh2')
  #now bfile and input should always be same sp just use bfile
  bfile=bfile+'_thresh2'
  #report on QC
  if (float(params['GENO']) < 1):
    proc = subprocess.Popen(['findstr', 'geno)', bfile+'.log'], stdout=subprocess.PIPE)
    tmp = proc.stdout.read()
    geno = int(tmp.split()[0])
  else:
    geno = 0
  logfile.write('Variants removed after --geno = '+str(geno)+'\n')
  proc = subprocess.Popen(['findstr', 'threshold(s)', bfile+'.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  maf = int(tmp.split()[0])

  logfile.write('Variants removed after --maf = '+str(maf)+'\n')
  proc = subprocess.Popen(['findstr', 'pass', bfile+'.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  oldv = v
  v = int(tmp.split()[0])
  logfile.write('\nvariants = '+str(v)+'\n')
  #check
  if oldv != (v+maf+geno):
    logfile.write('Check variant QC\n')
#2nd tier mind
  os.system(params['PLINKPATH']+' --bfile '+bfile+' --mind '+params['MIND']+' --make-bed --out '+bfile+'_thresh3')
  #now bfile and input should always be same sp just use bfile
  bfile=bfile+'_thresh3'
  #report on QC
  proc = subprocess.Popen(['findstr', 'mind)', bfile+'.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  mind = int(tmp.split()[0])
  logfile.write('Samples removed after --mind = '+str(mind)+'\n')
  proc = subprocess.Popen(['findstr', 'pass', bfile+'.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  oldn = n
  n = int(tmp.split()[3])
  logfile.write('\nn = '+str(n)+'\n')
  #check
  if oldn != (n+mind):
    logfile.write('Check sample QC\n')

#heterozygosity
#杂合性
if 'HET' in params:
  os.system(params['PLINKPATH']+' --bfile '+bfile+' --het --out '+bfile+'_het')
  os.system(params['RPATH']+' '+params['SCRIPTS']+'/het.R '+bfile+'_het.het > '+bfile+'_f.txt '+params['WORK'])
  #os.system('mv Rplots.pdf '+bfile+'_het.pdf')
  logfile.write('\nSamples sorted by f written to '+bfile+'_f.txt\n')
  hfile=open(bfile+'_hremove.txt', 'w')
  #if QC by F
  if 'FMIN' in params and 'FMAX' in params:
    fmin=float(params['FMIN'])
    fmax=float(params['FMAX'])
    fout=open(bfile+'_f.txt', 'r')
    first = 0
    hunder = 0
    hover = 0
    for line in fout:
      if first < 1:
        first = first + 1
      else:
        linesplit=line.split()
        if float(linesplit[6]) < fmin:
          hunder = hunder + 1
          hfile.write(linesplit[1]+'\t'+linesplit[2]+'\n')
        elif float(linesplit[6]) > fmax:
          hover = hover + 1
          hfile.write(linesplit[1]+'\t'+linesplit[2]+'\n')
    fout.close()

    logfile.write(str(hunder)+' sample(s) under FMIN\n'+str(hover)+' sample(s) over FMAX\n')
  #if QC by H
  else:
    hout=open(bfile+'_het.het', 'r')
    first = 0
    hunder = 0
    hover = 0
    hvector = []
    for line in hout:
      if first == 0:
        first = first + 1
      else:
        fid, iid, o_hom, e_hom, n_nm, f = line.split()
        hvector.append((float(n_nm)-float(o_hom))/float(n_nm))
    hout.close()
    hmean = statistics.mean(hvector)
    hstd = statistics.stdev(hvector)
    hout=open(bfile+'_het.het', 'r')
    first = 0
    for line in hout:
      if first == 0:
        first = first + 1
      else:
        fid, iid, o_hom, e_hom, n_nm, f = line.split()
        if ((float(n_nm)-float(o_hom))/float(n_nm)) < (hmean - 3*hstd):
          hunder = hunder +1
          hfile.write(fid+'\t'+iid+'\n')
        if ((float(n_nm)-float(o_hom))/float(n_nm)) > (hmean + 3*hstd):
          hover = hover +1
          hfile.write(fid+'\t'+iid+'\n')
    hout.close()
    logfile.write('H mean: '+str(hmean)+'\nH stdev: '+str(hstd)+'\n'+str(hunder)+' sample(s) under 3 stdev of H\n'+str(hover)+' sample(s) over 3 stdev of H\n')
  hfile.close()
  if hunder+hover > 0:
    os.system(params['PLINKPATH']+' --bfile '+bfile+' --remove '+bfile+'_hremove.txt --make-bed --out '+bfile+'_h')
    bfile=bfile+'_h'
    proc = subprocess.Popen(['findstr', 'remove:', bfile+'.log'], stdout=subprocess.PIPE)
    tmp = proc.stdout.read()
    n = int(tmp.split()[1])
    logfile.write('n =  '+str(n)+'\n\n')

#Hardy-Weinberg Equilibrium
#哈代-温伯格平衡
os.system(params['PLINKPATH']+' --bfile '+bfile+'  --hardy --pheno '+params['PFILE']+' --pheno-name '+params['PHENO']+' --out '+bfile+'_hwe')
os.system('findstr "O(HET)\|UNAFF" '+bfile+'_hwe.hwe > '+bfile+'_hwe_ctrl.hwe')
os.system(params['RPATH']+' '+params['SCRIPTS']+'/hwe.R '+bfile+'_hwe_ctrl.hwe '+params['WORK'])
#os.system('mv Rplots.pdf '+bfile+'_hwe.pdf')

hwe=float(params['HWE'])
hweout=open(bfile+'_hwe_ctrl.hwe', 'r')
hwefile=open(bfile+'_hweremove.txt', 'w')
first = 0
hweunder = 0
for line in hweout:
  if first < 1:
    first = first + 1
  else:
    linesplit=line.split()

    try:
      if float(linesplit[8]) < hwe:
        hweunder = hweunder + 1
        hwefile.write(linesplit[1]+'\n')
    except ValueError:
        logfile.write('Value Error in HWE: '+line)
hwefile.close()
hweout.close()
logfile.write('\nVariants with p-value lower than HWE: '+str(hweunder)+'\n')
if hweunder > 0:
  os.system(params['PLINKPATH']+' --bfile '+bfile+' --exclude '+bfile+'_hweremove.txt --make-bed --out '+bfile+'_hwe')
  bfile=bfile+'_hwe'
  proc = subprocess.Popen(['findstr', 'exclude:', bfile+'.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  v = int(tmp.split()[1])
  logfile.write('variants =  '+str(v)+'\n\n')

#结果打包到压缩文件
os.system('tar czvf '+params['WORK']+'/final_graphs.tgz '+params['WORK']+'/*.png')
logfile.write('Final bfile: '+bfile+'\nGraphs in final_graphs.tgz\n')

if 'CLEAN' in params:
  #os.system('mv '+bfile+'+.log final.log && mv '+bfile+'.bed final.bed.t && mv '+bfile+'.bim final.bim.t && mv '+bfile+'.fam final.fam.t && rm *.bed && rm *.bim && rm *.fam && mv final.bed.t final.bed && mv final.bim.t final.bim && mv final.fam.t final.fam ')
  os.system('mv '+bfile+'.log '+params['WORK']+'/final.log')
  os.system('mv '+bfile+'.bed '+params['WORK']+'/final.bed')
  os.system('mv '+bfile+'.bim '+params['WORK']+'/final.bim')
  os.system('mv '+bfile+'.fam '+params['WORK']+'/final.fam')
  os.system('rm '+params['INPUT']+'_*')
  print('\nFiles cleaned\n')

#判断n0的值是否为空
def check_step_N0(MAF,MIND,GENO,HWE):
  if MAF is not None:
    return 'Valid'
  else:
    return 'Invalid'
def check_step_N1(MAF,MIND,GENO,HWE):
  if MIND is not None:
    return 'Valid'
  else:
    return 'Invalid'
def check_step_N2(MAF,MIND,GENO,HWE):
  if GENO is not None:
    return 'Valid'
  else:
    return 'Invalid'
def check_step_N3(MAF,MIND,GENO,HWE):
  if HWE is not None:
    return 'Valid'
  else:
    return 'Invalid'

def smart_contract_validation(MAF,MIND,GENO,HWE):
  switcher = {0: check_step_N0(MAF,MIND,GENO,HWE),
              1: check_step_N1(MAF,MIND,GENO,HWE),
              2: check_step_N2(MAF,MIND,GENO,HWE),
              3: check_step_N3(MAF,MIND,GENO,HWE)
              }
  return switcher.get(MAF,"Valid")


post_status=post_blockchain_data(MAF,MIND,GENO,HWE,Individual,variants)
block_status = smart_contract_validation(MAF,MIND,GENO,HWE)
if block_status == 'Valid':
  print("QC is successful ")
else:
  print("QC is unsuccessful")
time_end = time.time()
print('time cost', time_end - time_start, 's')