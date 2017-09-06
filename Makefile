all:
	g++ --std=c++11 -pthread \
	src/main.cpp src/FastaHeuristicFilter.cpp src/NuclSeq.cpp \
	src/Alignment.cpp src/SeqProfile.cpp src/FileUtils.cpp \
	src/BigFasta.cpp \
	-o bin/alenhador
test:
	make all
	./bin/alenhador example-data/query.fa example-data/database.fa 3