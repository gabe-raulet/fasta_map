all: fasta_map.a

OBJS=mmapfile.o fasta_map.o

fasta_map.a: $(OBJS)
	ar rcs $@ $^

%.o: %.c
	gcc-11 -Wall -O2 -c -o $@ $<

clean:
	rm -rf $(OBJS)
