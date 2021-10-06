

@echo off

::========================
:: delete containers
::========================
set container_name_prefix=%1
set con_start_index=%2
set con_end_index=%3




(FOR /L %%G IN (%con_start_index%,1,%con_end_index%) DO (

	docker stop %container_name_prefix%_%%G
	docker container rm %container_name_prefix%_%%G


)
)