BUILD = build
SRC = src
INC = include
BIN = bin
INCLUDE = -I$(INC)
CC = g++ -g

CLASSES = $(BUILD)/ROC.o $(BUILD)/ScoreData.o $(BUILD)/getopt_pp.o $(BUILD)/Parameters.o

all : CVRLROC

CVRLROC : $(BUILD)/CVRLROC.o $(CLASSES)
	$(CC) -o $(BIN)/CVRLROC $(BUILD)/CVRLROC.o $(CLASSES) -static

$(BUILD)/CVRLROC.o : $(SRC)/CVRLROC.cpp
	$(CC) -c $(SRC)/CVRLROC.cpp $(INCLUDE)
	mv -f CVRLROC.o $(BUILD)/CVRLROC.o

$(BUILD)/AlgorithmComparison.o : $(SRC)/AlgorithmComparison.cpp
	$(CC) -c $(SRC)/AlgorithmComparison.cpp $(INCLUDE)
	mv -f AlgorithmComparison.o $(BUILD)/AlgorithmComparison.o

$(BUILD)/ROC.o : $(SRC)/ROC.cpp
	$(CC) -c $(SRC)/ROC.cpp $(INCLUDE)
	mv -f ROC.o $(BUILD)/ROC.o

$(BUILD)/ScoreData.o : $(SRC)/ScoreData.cpp
	$(CC) -c $(SRC)/ScoreData.cpp $(INCLUDE)
	mv -f ScoreData.o $(BUILD)/ScoreData.o

$(BUILD)/Parameters.o : $(SRC)/Parameters.cpp
	$(CC) -c $(SRC)/Parameters.cpp $(INCLUDE)
	mv -f Parameters.o $(BUILD)/Parameters.o

$(BUILD)/getopt_pp.o : $(SRC)/getopt_pp.cpp
	$(CC) -c $(SRC)/getopt_pp.cpp $(INCLUDE)
	mv -f getopt_pp.o $(BUILD)/getopt_pp.o

clean :
	rm -rf $(BUILD)/* $(BIN)/*
