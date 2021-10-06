@echo off

set container_name=%1
set index=%2

set data_file=%3
set cli_timeout=%4
set create_timeout=%5
set execution_timeout=%6
set control_level=%7


:: run individual container
docker start %container_name%
docker exec -it %container_name% /bin/bash /home/mythril/se_guider_run.sh %index% %data_file%  %cli_timeout% %create_timeout% %execution_timeout% %control_level%


:: stop container
docker stop %container_name%


:: delete container
docker container rm %container_name%

