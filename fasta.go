// Package fasta implements data structures and functions for reading, writing, and manipulating sequences in FASTA format.
package fasta

import (
	"bufio"
	"bytes"
	"io"
	"math"
	"math/rand"
)

const (
	DefaultLineLength = 70
)

// Sequence holds a nucleotide or protein sequence.
type Sequence struct {
	header     string
	data       []byte
	lineLength int
}

// A Sequence is read using a Scanner.
type Scanner struct {
	s                             *bufio.Scanner
	isHeader                      bool
	lastSequence                  bool
	previousHeader, currentHeader string
	firstSequence                 bool
	data                          []byte
}

func (s *Sequence) Header() string  { return s.header }
func (s *Sequence) Data() []byte    { return s.data }
func (s *Sequence) LineLength() int { return s.lineLength }

// If the line length passed is less than 1, it is assumed that effectively infinite lines are requested.
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
	f := []byte("ACGTUWSMKRYBDHVNacgtuwsmkrybdhvn")
	r := []byte("TGCAAWSKMYRVHDBNtgcaawskmyrvhdbn")
	var dic [256]byte
	for i, _ := range dic {
		dic[i] = byte(i)
	}
	for i, v := range f {
		dic[v] = r[i]
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

// ScanLine reads input line by line. It skips empty lines and marks headers.
func (s *Scanner) ScanLine() bool {
	for s.s.Scan() {
		b := s.s.Bytes()
		if len(b) > 0 {
			if b[0] == '>' {
				s.isHeader = true
			} else {
				s.isHeader = false
			}
			return true
		}
	}
	return false
}
func (s *Scanner) IsHeader() bool {
	return s.isHeader
}

// Line returns the last non-empty line scanned.
func (s *Scanner) Line() []byte {
	return s.s.Bytes()
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
	s.previousHeader = s.currentHeader
	if !s.firstSequence {
		return true
	} else {
		return false
	}
}

// NewScanner returns a new Scanner to read from r.
func NewScanner(r io.Reader) *Scanner {
	sc := bufio.NewScanner(r)
	scanner := Scanner{
		s:             sc,
		firstSequence: true,
	}
	return &scanner
}
