.PHONY: all bin clean

BIN = count_words
#DLL = counts.so

CFLAGS = -O2 -march=native

all: bin

bin: $(BIN)

count_words: count_words.o countscalc.o
	gcc $^ -o $@ -lz

count_words.o: count_words.c
	gcc $^ -c $(CFLAGS) -o $@

countscalc.o: countscalc.c
	gcc $^ -c $(CFLAGS) -o $@

#dll: $(DLL)
#
#counts.so: counts.o
#	gcc $^ -shared -o $@ -lz
#
#counts.o: countsmodule.c
#	gcc $^ -fPIC -c $(CFLAGS) -o $@

clean:
	rm -f *.pyc *.so *.o $(BIN)
