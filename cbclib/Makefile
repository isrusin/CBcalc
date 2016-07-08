.PHONY: all bin dll clean

SRC_DIR = src
BIN_DIR = bin
DLL_DIR = dll

BIN = $(BIN_DIR)/count_sw $(BIN_DIR)/count_lw
DLL = $(DLL_DIR)/count_words.so

CFLAGS = -O2 -march=native

all: bin dll

bin: $(BIN)

$(BIN_DIR)/count_%w: $(BIN_DIR)/count_%w.o $(BIN_DIR)/count_words.o
	gcc -o $@ -lz $^

$(BIN_DIR)/count_%.o: $(SRC_DIR)/count_%.c
	mkdir -p $(BIN_DIR)
	gcc -c $(CFLAGS) -o $@ $^

dll: $(DLL)

$(DLL_DIR)/count_words.so: $(DLL_DIR)/pic_count_words.o
	gcc -shared -o $@ -lz $^

$(DLL_DIR)/pic_count_words.o: $(SRC_DIR)/count_words.c
	mkdir -p $(DLL_DIR)
	gcc -fPIC -c $(CFLAGS) -o $@ $^

clean:
	rm -rf $(BIN_DIR) $(DLL_DIR)
	rm -f *.pyc
