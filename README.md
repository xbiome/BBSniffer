# BGCSniffer
## Running the pipeline
### Build your Docker image
在自己的计算机/集群构建镜像，命令参考下方，其中-t后引号中的内容为镜像名字和版本信息，用冒号分隔。
Build Docker image in users' PC or HPC cluster. Command is shown blow, content quoted by single quote after '-t' is the name and version info of Docker image, separated by ":".
```
docker build -t 'bgcsniffer:1' .
```
### run program in command line
镜像构建完成后运行整个流程，命令见下：
Run the pipeline with following command:
```
docker run -i --rm -v WORK_DIR:/in  bgcsniffer:1 bash -c "python /opt/DockerImage/main.py -json /in/test.json -workdir /in/"
```
用户需将上方命令中的WORK_DIR替换为自己的分析路径
Users should substitute 'WORK_DIR' with their own working directory, replace 'bgcsniffer:1' with their defined info from previous step and change the parameters in test.json file. Other command should be remained the same.
