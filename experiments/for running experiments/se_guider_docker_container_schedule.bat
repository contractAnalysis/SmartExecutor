

@echo off

::========================
:: create containers
::========================
set container_name_prefix=%1
set con_start_index=%2
set con_end_index=%3

set cli_timeout=%4
set create_timeout=%5
set execution_timeout=%6

set control_level=%7

set image=1b615fb0bf9d
set workspace=C:\Users\SERC\wei_experiments\experiment_paper\
set run_one_container=se_guider_run_for_one_container_window.bat

:: pay attention to the prefix of data file
set data_file_prefix=contracts_info_part

::========================
:: create containers
::========================

(FOR /L %%G IN (%con_start_index%,1,%con_end_index%) DO (

	echo create container: %container_name_prefix%_%%G	
	docker run -it --cpus 2 -m 16384m -v %workspace%:/home/mythril --ulimit stack=100000000:100000000 --name %container_name_prefix%_%%G --entrypoint /bin/bash %image% 

	
)
)


::========================
:: start docker ontainers
::========================
(FOR /L %%G IN (%con_start_index%,1,%con_end_index%) DO (

	:: arguments order: container, index, data file, cli_timeout, create_timeout,execution_timeout
 	echo %container_name_prefix%_%%G %%G %data_file_prefix%_%%G_r.csv %cli_timeout% %create_timeout% %execution_timeout% %control_level%
	
	:: pay attention to *.csv file: have suffix '_r' or not 
	start %run_one_container% %container_name_prefix%_%%G %%G %data_file_prefix%_%%G_r.csv %cli_timeout% %create_timeout% %execution_timeout% %control_level%

	

)
)
	
	

