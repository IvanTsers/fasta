// Package fasta implements data structures and functions for reading, writing, and manipulating sequences in FASTA format.
package fasta

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"math"
	"math/rand"
	"os"
)

const (
	DefaultLineLength = 70
)

var dic []byte

// Sequence holds a nucleotide or protein sequence.
type Sequence struct {
	header     string
	data       []byte
	lineLength int
}

// A Sequence is read using a Scanner.
type Scanner struct {
	r                             *bufio.Reader
	line                          []byte
	err                           error
	isHeader                      bool
	lastSequence                  bool
	previousHeader, currentHeader string
	firstSequence                 bool
	data                          []byte
}

func (s *Sequence) Header() string  { return s.header }
func (s *Sequence) Data() []byte    { return s.data }
func (s *Sequence) LineLength() int { return s.lineLength }

// SetHeader replaces the existing header.
func (s *Sequence) SetHeader(h string) {
	s.header = h
}

// SetData replaces the existing data.
func (s *Sequence) SetData(d []byte) {
	s.data = d
}

// SetLineLength replaces the current line length. If the line length passed is less than 1, it is assumed that effectively infinite lines are requested.
func (s *Sequence) SetLineLength(l int) {
	s.lineLength = l
	if s.lineLength < 1 {
		s.lineLength = math.MaxInt64
	}
}

// AppendToHeader appends the suffix suf to the header.
func (s *Sequence) AppendToHeader(suf string) {
	s.header = s.header + suf
}

// Equals compares two sequences and returns true if their headers and data are identical.
func (a *Sequence) Equals(b *Sequence) bool {
	if a.header != b.header {
		return false
	}
	if len(a.data) != len(b.data) {
		return false
	}
	return bytes.Equal(a.data, b.data)
	return true
}

// String wraps the sequence into lines at most lineLength characters long.
func (s *Sequence) String() string {
	var b []byte
	b = append(b, '>')
	b = append(b, s.header...)
	b = append(b, '\n')
	var c int
	for _, r := range s.data {
		b = append(b, r)
		c++
		if c == s.lineLength {
			c = 0
			b = append(b, '\n')
		}
	}
	if c == 0 && len(b) > 0 {
		b = b[:len(b)-1]
	}
	if len(b) > 0 {
		return string(b)
	} else {
		return ""
	}
}

// Method Shuffle randomizes the residues in a Sequence. The sequence composition remains unchanged.
func (s *Sequence) Shuffle(r *rand.Rand) {
	d := s.data
	r.Shuffle(len(d), func(i, j int) {
		d[i], d[j] = d[j], d[i]
	})
}

// Method Reverse reverses the residues of a Sequence.
func (s *Sequence) Reverse() {
	d := s.data
	for i, j := 0, len(d)-1; i < j; i, j = i+1, j-1 {
		d[i], d[j] = d[j], d[i]
	}
}

// Complement complements nucleotide sequences.
func (s *Sequence) Complement() {
	if dic == nil {
		dic = make([]byte, 256)
		f := []byte("ACGTUWSMKRYBDHVNacgtuwsmkrybdhvn")
		r := []byte("TGCAAWSKMYRVHDBNtgcaawskmyrvhdbn")
		for i, _ := range dic {
			dic[i] = byte(i)
		}
		for i, v := range f {
			dic[v] = r[i]
		}
	}
	for i, v := range s.data {
		s.data[i] = dic[v]
	}
}

// ReverseComplement reverse-complements a Sequence.
func (s *Sequence) ReverseComplement() {
	s.Reverse()
	s.Complement()
}

// Method Length returns the number of residues in Sequence.
func (s *Sequence) Length() int {
	return len(s.data)
}

// Method GC returns the fraction of GC nucleotides in Sequence.
func (s *Sequence) GC() float64 {
	l := float64(s.Length())
	gc := 0.0
	for _, r := range s.data {
		if r == 'G' || r == 'C' {
			gc++
		}
	}
	return gc / l
}

// Method Clean removes non-canonical nucleotides from a Sequence (that is, keeps only ATGC/atgc).
func (s *Sequence) Clean() {
	d := s.Data()
	i := 0
	for _, c := range d {
		if c == 'A' || c == 'C' || c == 'G' || c == 'T' ||
			c == 'a' || c == 'c' || c == 'g' || c == 't' {
			d[i] = c
			i++
		}
	}
	d = d[:i]
	s.SetData(d)
}

// Method DataToUpper converts Data bytes to uppercase.
func (s *Sequence) DataToUpper() {
	d := s.Data()
	d = bytes.ToUpper(d)
	s.SetData(d)
}

// ScanLine reads input line by line. It skips empty lines and marks headers. The last call to ScanLine should be followed by a call to Flush to retrieve any bytes not terminated by newline.
func (s *Scanner) ScanLine() bool {
	var err error
	s.line, err = s.r.ReadBytes('\n')
	if err != nil {
		s.err = err
		return false
	}
	s.line = bytes.TrimRight(s.line, "\r\n")
	if len(s.line) > 0 {
		if s.line[0] == '>' {
			s.isHeader = true
		} else {
			s.isHeader = false
		}
		return true
	}
	s.err = nil
	return true
}
func (s *Scanner) IsHeader() bool {
	return s.isHeader
}

// Line returns the last non-empty line scanned.
func (s *Scanner) Line() []byte {
	return s.line
}

//  Flush returns any bytes remaining in the buffer after the  last call to ScanLine.
func (s *Scanner) Flush() []byte {
	var dum []byte
	if s.err == io.EOF {
		return s.line
	}
	return dum
}

// Sequence returns the last Sequence scanned.
func (s *Scanner) Sequence() *Sequence {
	seq := &Sequence{
		header: s.previousHeader,
	}
	seq.data = make([]byte, len(s.data))
	copy(seq.data, s.data)
	seq.lineLength = DefaultLineLength
	s.data = s.data[:0]
	return seq
}

// Function NewSequence returns a new Sequence.
func NewSequence(h string, d []byte) *Sequence {
	s := new(Sequence)
	s.header = h
	s.data = make([]byte, len(d))
	copy(s.data, d)
	s.lineLength = DefaultLineLength
	return s
}

//  ScanSequence reads input Sequence by Sequence.
func (s *Scanner) ScanSequence() bool {
	if s.lastSequence {
		return false
	}
	for s.ScanLine() {
		if s.isHeader {
			s.previousHeader = s.currentHeader
			s.currentHeader = string(s.Line()[1:])
			if s.firstSequence {
				s.firstSequence = false
			} else {
				return true
			}
		} else {
			s.data = append(s.data, s.Line()...)
		}
	}
	s.lastSequence = true
	if s.err == io.EOF {
		s.data = append(s.data, s.Line()...)
	}
	s.previousHeader = s.currentHeader
	if !s.firstSequence {
		return true
	} else {
		return false
	}
}

// NewScanner returns a new Scanner to read from r.
func NewScanner(r io.Reader) *Scanner {
	rd := bufio.NewReader(r)
	scanner := Scanner{
		r:             rd,
		firstSequence: true,
	}
	return &scanner
}

// ReadAll reads all sequences from a file and returns a slice of Sequences.
func ReadAll(f *os.File) []*Sequence {
	sc := NewScanner(f)
	var s []*Sequence
	for sc.ScanSequence() {
		s = append(s, sc.Sequence())
	}
	f.Close()
	return s
}

// Concatenate accepts a slice of Sequences and a sentinel byte. It concatenates the slice into a single Sequence entry, where all headers and data are glued. The concatenated headers and pieces of data are separated with the sentinel byte, if the latter is not zero.
func Concatenate(seqSlice []*Sequence, sentinel byte) *Sequence {
	l := len(seqSlice)
	switch {
	case l > 1:
		h := []byte(seqSlice[0].Header())
		d := seqSlice[0].Data()
		for i := 1; i < l; i++ {
			if sentinel != 0 {
				h = append(h, sentinel)
				d = append(d, sentinel)
			}
			h = append(h, []byte(seqSlice[i].Header())...)
			d = append(d, seqSlice[i].Data()...)
		}
		cSeq := NewSequence(string(h), d)
		return cSeq
	case l == 1:
		return seqSlice[0]
	default:
		fmt.Fprintln(os.Stderr,
			"fasta.Concatenate: the input slice is empty")
		return nil
	}
}
