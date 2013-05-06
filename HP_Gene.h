#ifndef _HP_GENE
#define _HP_GENE

#include <list>
#include <string>
#include <fstream>

using namespace std;

class HP_Record {
public:
	/*
	The ID of the landmark used to establish the coordinate system for the current feature. IDs may contain any characters, but must escape any characters not in the set [a-zA-Z0-9.:^*$@!+_?-|]. In particular, IDs may not contain unescaped whitespace and must not begin with an unescaped ">".
	To escape a character in this, or any of the other GFF3 fields, replace it with the percent sign followed by its hexadecimal representation. For example, ">" becomes "%E3".	   
	*/
	string seqid;
	/*
	The source is a free text qualifier intended to describe the algorithm or operating procedure that generated this feature. Typically this is the name of a piece of software, such as "Genescan" or a database name, such as "Genbank." In effect, the source is used to extend the feature ontology by adding a qualifier to the type creating a new composite type that is a subclass of the type in the type column. It is not necessary to specify a source. If there is no source, put a "." (a period) in this field.
	*/
	string source;
	/*
	The type of the feature (previously called the "method"). This is constrained to be either: (a) a term from the "lite" sequence ontology, SOFA; or (b) a SOFA accession number. The latter alternative is distinguished using the syntax SO:000000. This field is required.
	*/
	string type;
	/*
	The start and end of the feature, in 1-based integer coordinates, relative to the landmark given in column 1. Start is always less than or equal to end.
	For zero-length features, such as insertion sites, start equals end and the implied site is to the right of the indicated base in the direction of the landmark. These fields are required.
	*/
	int start;
	int end;
	/*
	The score of the feature, a floating point number. As in earlier versions of the format, the semantics of the score are ill-defined. It is strongly recommended that E-values be used for sequence similarity features, and that P-values be used for ab initio gene prediction features. If there is no score, put a "." (a period) in this field.   
	*/
	float score;
	/*
	The strand of the feature. + for positive strand (relative to the landmark), - for minus strand, and . for features that are not stranded. In addition, ? can be used for features whose strandedness is relevant, but unknown.
	*/
	char strand;
	/*
	For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame. The phase is one of the integers 0, 1, or 2, indicating the number of bases that should be removed from the beginning of this feature to reach the first base of the next codon. In other words, a phase of "0" indicates that the next codon begins at the first base of the region described by the current line, a phase of "1" indicates that the next codon begins at the second base of this region, and a phase of "2" indicates that the codon begins at the third base of this region. This is NOT to be confused with the frame, which is simply start modulo 3.

	For forward strand features, phase is counted from the start field. For reverse strand features, phase is counted from the end field.

	The phase is REQUIRED for all CDS features.
	*/
	int phase;
	/*
	Indicates the unique identifier of the feature. IDs must be unique within the scope of the GFF file.
	*/
	string ID;
	/*
	Display name for the feature. This is the name to be displayed to the user. Unlike IDs, there is no requirement that the Name be unique within the file.
	*/
	string name;
	/*
	Indicates the parent of the feature. A parent ID can be used to group exons into transcripts, transcripts into genes, an so forth. A feature may have multiple parents. Parent can *only* be used to indicate a partof relationship.
	*/
	string parent;

	// constructor
	HP_Record();
	string toString() const;
};

class HP_CDS: public HP_Record {
};

class HP_Exon: public HP_Record {
public:
	//list<CDS> cdss;
	string toString() const;
};

class HP_MRNA: public HP_Record {
public:
	list<HP_Exon> exons;
	//list<CDS> cdss;
	string toString() const;
	int getLength() const;
};

class HP_Gene: public HP_Record {
public:
	list<HP_MRNA> mRNAs;
	string toString() const;
	void getBounds(int &beg, int&end) const;
};

#endif
