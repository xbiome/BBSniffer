# BGCSniffer
## Build your Docker image
```
docker build -t 'bgcsniffer:1' .
```
run program
```
docker run -i --rm -v WORK_DIR:/in  bgcsniffer:1 bash -c "python /opt/DockerImage/main.py -json /in/test.json -workdir /in/"
```
