# BGCSniffer
## Running the pipeline

### Creating Input file
根据项目需求更改输入文件test.json。输入文件参考格式请见下:<br>
Edit input file base on your own project. Input file format was provided as follow:
```
{
    "Query": "(pf16570 OR pf16569 OR pf16555) and (ipr042002)",
    "Input_genomes": "Input_genomes",
    "Email": "zhuzhengnong@xbiome.com",
    "Antismash_threads": 10,
    "Antismash_gap_len": 20,
    "Antismash_extenson_len": 8
}
```
1. Query<br>
	Input Protein family ID. Same category of protein families should be included in one parenthesis and linked by ' OR ', Different categories of protein families shold be linked with ' and '.
2. Input_genomes<br>
	Name of folder which contains input reference model strain's complete genome file(s).
3. Email<br>
	Email of users.
4. Antismash_threads<br>
	Thread number when runnning antismash. Higher thread number may reduce run time but may cost more computing resources.
5. Antismash_gap_len<br>
	Maximum length(in kbp) of gap region between elements in one BGC.
6. Antismash_extenson_len<br>
	Length(in kbp) of flanking region of BGC's core region.

### Build your Docker image
在自己的计算机/集群构建镜像，命令参考下方，其中-t后引号中的内容为镜像名字和版本信息，用冒号分隔。<br>
Build Docker image in users' PC or HPC cluster. Command is shown blow, content quoted by single quote after '-t' is the name and version info of Docker image, separated by ":".
```
docker build -t 'bmpsniffer:1' .
```
### Run program in command line
镜像构建完成后运行整个流程，命令见下：<br>
Run the pipeline with following command:
```
docker run -i --rm -v WORK_DIR:/in  bmpsniffer:1 bash -c "python /opt/DockerImage/main.py -json /in/test.json -workdir /in/"
```
用户需将上方命令中的WORK_DIR替换为自己的分析路径。<br>
Users should substitute 'WORK_DIR' with their own working directory, replace 'bmpsniffer:1' with their defined info from previous step and change the parameters in test.json file. Other command should be remained the same.

### Run Protein Keywords searching
用户可用此脚本进行蛋白质库中的Keywords searching。<br>
This scritp could be used to search Keywords in protein database (PFAM) and acquire relevant protein family ID.
```
docker run -i --rm -v WORK_DIR:/in  bmpsniffer:1 bash -c "python /opt/DockerImage/KeywordsSearch.py -work_dir /in/ -input_keywords 'gram positive pilin'"
```
运行结束后可在工作目录下查看ProteinID.list文件。
