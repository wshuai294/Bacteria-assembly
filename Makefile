all:	select_ref	continuous
select_ref: scripts/src/find_best_assembly.cpp
	g++ -pthread -o scripts/src/select_ref scripts/src/find_best_assembly.cpp
continuous: scripts/src/find_continuous_seg.cpp
	g++ -pthread -o scripts/src/continuous scripts/src/find_continuous_seg.cpp
