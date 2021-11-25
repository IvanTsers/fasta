package fasta

import (
	"bytes"
	"fmt"
	"io/ioutil"
	"math/rand"
	"os"
	"strconv"
	"testing"
)

func TestEquals(t *testing.T) {
	d1 := []byte("ACCGT")
	s1 := NewSequence("s1", d1)
	s2 := NewSequence("s1", d1)
	if !(s1.Equals(s2) && s2.Equals(s1)) {
		t.Error("equal sequences declared unequal")
	}
	s3 := NewSequence("s3", d1)
	if s1.Equals(s3) || s3.Equals(s1) {
		t.Error("sequences with unequal headers declared equal")
	}
	d2 := []byte("ACGGT")
	s4 := NewSequence("s1", d2)
	if s1.Equals(s4) || s4.Equals(s1) {
		t.Error("sequences with unequal data declared equal")
	}
}
func TestString(t *testing.T) {
	seq := NewSequence("seq", []byte("ACCGT"))
	ll := []int{4, 5, 10, 0}
	for _, l := range ll {
		seq.SetLineLength(l)
		f, _ := ioutil.TempFile("", "test_*")
		fmt.Fprintf(f, "%s\n", seq)
		f.Close()
		b1, _ := ioutil.ReadFile(f.Name())
		b2 := bytes.Split(b1, []byte("\n"))
		os.Remove(f.Name())
		h := string(b2[0][1:])
		if h != seq.Header() {
			t.Errorf("did not write header correctly; "+
				"want %q; get %q\n", h, seq.Header())
		}
		var counter int
		for i, s := range b2 {
			if i < 1 {
				continue
			}
			for _, c := range s {
				d := seq.Data()
				if c != d[counter] {
					t.Error("did not write data correctly")
				}
				counter++
			}
		}
		if l <= 5 && l > 0 {
			if len(b2[1]) != l {
				t.Errorf("did not format data correctly: "+
					"want: %d; get: %d", l, len(b2[1]))
			}
		}
	}
}
func TestShuffle(t *testing.T) {
	orig := NewSequence("", []byte("ACCGT"))
	shuf := []byte("GTACC")
	r := rand.New(rand.NewSource(13))
	orig.Shuffle(r)
	if !bytes.Equal(orig.data, shuf) {
		t.Errorf("want:\n%s\nget:\n%s\n",
			string(shuf),
			string(orig.data))
	}
}
func TestReverse(t *testing.T) {
	ori := NewSequence("", []byte("ACCGT"))
	rev := []byte("TGCCA")
	ori.Reverse()
	if !bytes.Equal(ori.data, rev) {
		t.Errorf("want:\n%s\nget:\n%s\n",
			string(rev),
			string(ori.data))
	}
}
func TestReverseComplement(t *testing.T) {
	ori := NewSequence("", []byte("ACCGT"))
	rc := []byte("ACGGT")
	ori.ReverseComplement()
	if !bytes.Equal(ori.data, rc) {
		t.Errorf("want:\n%s\nget:\n%s\n",
			string(rc),
			string(ori.data))
	}
}

func TestScanner(t *testing.T) {
	for i := 1; i <= 8; i++ {
		name := "./data/seq" + strconv.Itoa(i) + ".fasta"
		in, err := os.Open(name)
		if err != nil {
			t.Errorf("couldn't open %q\n", name)
		}
		out, _ := ioutil.TempFile(".", "test_*")
		scanner := NewScanner(in)
		var foundSequence = false
		for scanner.ScanSequence() {
			foundSequence = true
			seq := scanner.Sequence()
			fmt.Fprintf(out, "%s\n", seq)
		}
		in.Close()
		out.Close()
		if foundSequence {
			id, _ := ioutil.ReadFile(name)
			od, _ := ioutil.ReadFile(out.Name())
			if !bytes.Equal(id, od) {
				t.Errorf("failed to reproduce %q\n", name)
			}
		}
		os.Remove(out.Name())
	}
}
