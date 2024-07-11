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
func TestLength(t *testing.T) {
	nuc := "ACCGT"
	seq := NewSequence("", []byte(nuc))
	l := seq.Length()
	if l != len(nuc) {
		t.Errorf("want:\n%d\nget:\n%d\n",
			len(nuc), l)
	}
}
func TestGC(t *testing.T) {
	s := []string{"ACCGT", "GGC", "AATAT"}
	w := []float64{3.0 / 5.0, 3.0 / 3.0, 0.0 / 5.0}
	g := 1.1
	for i, r := range s {
		seq := NewSequence("", []byte(r))
		g = seq.GC()
		if g != w[i] {
			t.Errorf("want:\n%v\nget:\n%v\n",
				w[i], g)
		}
	}
}
func TestKeepOnlyATGC(t *testing.T) {
	s := "XXATATNGTnCactAploenTTg"
	w := "ATATGTCactATTg"
	seq := NewSequence("", []byte(s))
	seq.KeepOnlyATGC()
	g := string(seq.Data())
	if g != w {
		t.Errorf("want:\n%s\nget:\n%s\n", w, g)
	}
}
func TestScanner(t *testing.T) {
	for i := 1; i <= 9; i++ {
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
			if i == 9 {
				id = append(id, '\n')
			}
			od, _ := ioutil.ReadFile(out.Name())
			if !bytes.Equal(id, od) {
				t.Errorf("failed to reproduce %q\n", name)
			}
		}
		os.Remove(out.Name())
	}
}
func TestFlush(t *testing.T) {
	f, _ := os.Open("data/seq8.fasta")
	sc := NewScanner(f)
	w := 5085
	g := 0
	for sc.ScanLine() {
		g += len(sc.Line())
	}
	g += len(sc.Flush())
	if g != w {
		t.Errorf("get:\n%d\nwant:\n%d\n", g, w)
	}
	f.Close()
	f, _ = os.Open("data/seq9.fasta")
	sc = NewScanner(f)
	g = 0
	for sc.ScanLine() {
		g += len(sc.Line())
	}
	g += len(sc.Flush())
	if g != w {
		t.Errorf("get:\n%d\nwant:\n%d\n", g, w)
	}
	f.Close()
}
func TestReadAndConcatenate(t *testing.T) {
	wantDataLen := [9]int{0, 0, 0, 5, 2 * 70, 2 * 140, 5 * 700, 5 * 1000,
		5 * 1000}
	wantHeaders := [9]string{"",
		"|",
		"|seq3",
		"|Rand_1; G/C=0.20",
		"|Rand_1; G/C=0.41|Rand_2; G/C=0.54",
		"|Rand_1; G/C=0.50|Rand_2; G/C=0.50",
		"|Rand_1; G/C=0.50|Rand_2; G/C=0.50|Rand_3; " +
			"G/C=0.50|Rand_4; G/C=0.50|Rand_5; G/C=0.50",
		"|Rand_1; G/C=0.50|Rand_2; G/C=0.50|Rand_3; " +
			"G/C=0.50|Rand_4; G/C=0.50|Rand_5; G/C=0.50",
		"|Rand_1; G/C=0.50|Rand_2; G/C=0.50|Rand_3; " +
			"G/C=0.50|Rand_4; G/C=0.50|Rand_5; G/C=0.50"}
	for i := 1; i <= 9; i++ {
		name := "./data/seq" + strconv.Itoa(i) +
			".fasta"
		f, err := os.Open(name)
		if err != nil {
			t.Errorf("couldn't open %q\n", name)
		}
		seq := ReadAndConcatenate(f)
		wl := wantDataLen[i-1]
		gl := seq.Length()
		if gl != wl {
			t.Errorf("Data:\nget:\n%d\nwant:\n%d\n", gl, wl)
		}
		wh := wantHeaders[i-1]
		gh := seq.Header()
		if gh != wh {
			t.Errorf("Headers:\nget:\n%s\nwant:\n%s\n", gh, wh)
		}
	}
}
