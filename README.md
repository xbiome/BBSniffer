# BBSniffer

BBSniffer is a software designed to identify suitable chassis strains from a comprehensive strain library, which encompasses pathogens, industrial microorganisms, and various non-pathogenic species, using specified biopolymer terms. It accomplishes this by analyzing biopolymer rules, HMM profiles, and genomic sequences derived from the user's input.
![Description](https://github.com/xbiome/BBSniffer/blob/develop/Description.jpeg)


## 1 Data Download and Extraction
Use the following DOI links to download the data files from Zenodo (https://zenodo.org/), and untar the downloaded data files to a directory

   - 10.5281/zenodo.8286952
   - 10.5281/zenodo.8313016
   - 10.5281/zenodo.8318303
   - 10.5281/zenodo.8321312
   - 10.5281/zenodo.8321340
   - 10.5281/zenodo.8325104

## 2 Running the pipeline

### 2.1 Creating Input file
Edit input file according to your project. Input file format was provided as below:

```
{
    "Query": "(pf16570 OR pf16569 OR pf16555) and (ipr042002)",
    "Input_genomes": "Input_genomes",
    "Email": "somebody@123.com",
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

### 2.2 Clone git repository
Clone git repository of BBSniffer.

```
git clone https://github.com/xbiome/BBSniffer.git
cd BBSniffer
```
###

### 2.3 Build your Docker image
Build Docker image in users' PC or HPC cluster. Command is shown blow, content quoted by single quote after '-t' is the name and version info of Docker image, separated by ":".

```
docker build -t 'bbsniffer:1' .
```
### 2.4 Run program in command line
Run the pipeline with following command:

```
docker run -i --rm -v WORK_DIR:/in  bbsniffer:1 bash -c "python /opt/DockerImage/main.py -json /in/test.json -workdir /in/ --database /path/to/your/database/"
```
Users should substitute 'WORK_DIR' with their own working directory, replace 'bbsniffer:1' with their defined info from previous step and change the parameters in test.json file. Other command should be remained the same.

### 2.5 Run Protein Keywords searching
This scritp could be used to search Keywords in protein database (PFAM) and acquire relevant protein family ID.

```
docker run -i --rm -v WORK_DIR:/in  bbsniffer:1 bash -c "python /opt/DockerImage/KeywordsSearch.py -work_dir /in/ -input_keywords 'gram positive pilin'"
```
Check ProteinID.list file in working directory after the finish of the pipeline.

## 3. Tutorial for Windows users

### 3.1 Installation of WSL2
For Windows Users, the installation of WSL(Windows Subsystem for Linux) is required. <br>
Follow [this link](https://docs.microsoft.com/zh-cn/windows/wsl/install-manual) for complete installation of WSL(Windows Subsystem for Linux) on Windows 10 system PC. <br>
:+1: **Note: Ubuntu is recommanded in Step6** :+1:

### 3.2 Installation of latest Docker Desktop
Follow [this link](https://docs.docker.com/desktop/windows/install/) for complete installation of Docker Desktop.

### 3.3 Running BBSniffer on WSL
Open Ubuntu from your software list, login and Run pipeline following the instruction in Session 1


## 4. System requirement
### 4.1 Linux
- Avaiable space of system hard drive should be no less than 200G.
### 4.2 Windows
- 1903 or higher version are required for x64 system.
- 2004 or higher version are required for ARM system.
- WSL 2 is not supported on system lower than 18362.
- Avaiable space of system hard drive should be no less than 200G.
