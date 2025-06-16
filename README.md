# HORSCAN

**HORSCAN** is an innovative algorithm based on multi-level dynamic programming, specifically designed to identify structural variants at the Higher-Order Repeat (HOR) level in human centromeric regions.

---

## ⚙️ Installation

### Prerequisites

Ensure the [Rust toolchain](https://rustup.rs/) (including `rustc` and `cargo`) is installed on your system.

### Build from Source

Clone the repository and build the project with a single set of commands:

```bash
git clone [https://github.com/XDwan/HORSCAN.git](https://github.com/XDwan/HORSCAN.git)
cd HORSCAN
cargo build --release
```

The optimized executable will be available at `./target/release/HORSCAN`.

For system-wide access, you can move the executable to a directory in your `$PATH`:

```bash
sudo mv ./target/release/HORSCAN /usr/local/bin/
```

## 🚀 Usage

### Command-line Arguments

| Option                          | Shorthand                   | Description                                                      |
| :------------------------------ | :-------------------------- | :--------------------------------------------------------------- |
| `--source <SOURCE>`           | `-s <SOURCE>`             | **[Required]**Path to the source input BED file.                 |
| `--target <TARGET>`           | `-t <TARGET>`             | **[Required]**Path to the target input BED file.                 |
| `--output <OUTPUT>`           | `-o <OUTPUT>`             | **[Required]**Path for the output alignment file.                |
| `--mode <MATCH MISMATCH GAP>` | `-m <MATCH MISMATCH GAP>` | **[Required]**Set scores for Match, Mismatch, and Gap penalties. |

### 📄 Input File Format

Input files must be **tab-separated** text files in a 4-column BED-like format.

1. **Sample/Chromosome** : Identifier for the sequence (e.g., `SampleA#chrX`).
2. **Start** : Start coordinate (integer, 0-based).
3. **End** : End coordinate (integer).
4. **Monomer Label** : A string identifying the monomer unit.

### Full Example

1. **Create example input files:**

   ```bash
   # Create source.bed (Sequence: A B D E F G H I J K)
   echo -e "CHM13#chrX\t1000\t1170\tA\nCHM13#chrX\t1171\t1341\tB\nCHM13#chrX\t1342\t1512\tD\nCHM13#chrX\t1513\t1683\tE\nCHM13#chrX\t1684\t1854\tF\nCHM13#chrX\t1855\t2025\tG\nCHM13#chrX\t2026\t2196\tH\nCHM13#chrX\t2197\t2367\tI\nCHM13#chrX\t2368\t2538\tJ\nCHM13#chrX\t2539\t2709\tK" > source.bed

   # Create target.bed (Sequence: A B C D E O G H J K)
   echo -e "CHM1#chrX\t1000\t1170\tA\nCHM1#chrX\t1171\t1341\tB\nCHM1#chrX\t1342\t1512\tC\nCHM1#chrX\t1513\t1683\tD\nCHM1#chrX\t1684\t1854\tE\nCHM1#chrX\t1855\t2025\tO\nCHM1#chrX\t2026\t2196\tG\nCHM1#chrX\t2197\t2367\tH\nCHM1#chrX\t2368\t2538\tJ\nCHM1#chrX\t2539\t2709\tK" > target.bed
   ```
2. **Run HORSCAN:**

   ```bash
   # If HORSCAN is in your PATH
   HORSCAN --source source.bed --target target.bed --output test --mode 10 4 2

   # If running directly from the project directory
   ./target/release/HORSCAN -s source.bed -t target.bed -o test -m 10 4 2
   ```

### 📋 Output File Format

The output is a **tab-separated** file describing the pairwise alignment of monomers, with 9 columns.

| Columns | Content          | Description                                                 |
| :------ | :--------------- | :---------------------------------------------------------- |
| 1-4     | Source Monomer   | `Sample`,`Start`,`End`,`Label`from the source file. |
| 5-8     | Target Monomer   | `Sample`,`Start`,`End`,`Label`from the target file. |
| 9       | Alignment Status | A code:`MTH`,`MIS`,`INS`,`DEL`,.                    |



### Example Output Explained

Running the command above will produce an output file `test.alignment` with content similar to this:

```
CHM13#chrX	1000	1170	A	CHM1#chrX	1000	1170	A	MTH
CHM13#chrX	1171	1341	B	CHM1#chrX	1171	1341	B	MTH
CHM13#chrX	1171	1171	-	CHM1#chrX	1342	1512	C	DEL
CHM13#chrX	1342	1512	D	CHM1#chrX	1513	1683	D	MTH
CHM13#chrX	1513	1683	E	CHM1#chrX	1684	1854	E	MTH
CHM13#chrX	1684	1854	F	CHM1#chrX	1855	2025	O	MIS
CHM13#chrX	1855	2025	G	CHM1#chrX	2026	2196	G	MTH
CHM13#chrX	2026	2196	H	CHM1#chrX	2197	2367	H	MTH
CHM13#chrX	2197	2367	I	CHM1#chrX	2197	2197	-	INS
CHM13#chrX	2368	2538	J	CHM1#chrX	2368	2538	J	MTH
CHM13#chrX	2539	2709	K	CHM1#chrX	2539	2709	K	MTH
```

* **Line 3 (Insertion):** A gap in the source is aligned to monomer `C` from the target, indicating an insertion of `C` in the target sequence.
* **Line 6 (Mismatch):** Monomer `F` in the source aligns to a different monomer `O` in the target.
* **Line 9 (Deletion):** Monomer `I` in the source is aligned to a gap in the target, indicating a deletion of `I` from the target sequence.

---

## 💡 Future Improvements

We are actively working on enhancing HORSCAN. Planned features include:

* **Visualization Module** : To generate dot plots and linear diagrams for intuitive representation of HOR alignments and variations.
* **Seamless Integration** : Direct support for parsing outputs from standard centromere annotation tools like HiCAT and HORmon.
* **VCF Output Format** : Implementation of VCF output for HOR-level structural variants, improving compatibility with downstream analysis pipelines.
* **Performance Optimization** : Further performance gains through multi-threading and algorithmic enhancements.
* **Enhanced HOR-aware Alignment** : Refinement of the core algorithm to better handle complex variations and improve annotation consistency across alignments.
