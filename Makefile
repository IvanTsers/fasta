all : fasta

fasta: fasta.go
	go build fasta.go
fasta.go: fasta.org
	awk -f scripts/preTangle.awk fasta.org | bash scripts/org2nw | notangle -Rfasta.go | gofmt > fasta.go
test: fasta_test.go fasta.go
	go test -v
fasta_test.go: fasta_test.org
	awk -f scripts/preTangle.awk fasta_test.org | bash scripts/org2nw | notangle -Rfasta_test.go | gofmt > fasta_test.go

.PHONY: doc
doc:
	make -C doc

clean:
	rm -f *.go
	make clean -C doc
