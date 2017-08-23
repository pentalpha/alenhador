all:
	g++ --std=c++11 src/main.cpp src/NuclSeq.cpp src/Alignment.cpp src/SeqProfile.cpp src/FileUtils.cpp -o bin/alenhador
test:
	./bin/alenhador example-data/query.fa example-data/database.fa 6