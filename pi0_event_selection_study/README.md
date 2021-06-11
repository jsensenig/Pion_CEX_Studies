### Charge Exchange Event Selection.

This assumes at least Root v6.22 and Cmake 3.15 but should work with slightly older versions, just make the change in 
the `CMakeLists.txt` file. To compile do the following,
```python
mkdir build
```
and set up using cmake
```python
cmake . -B build
```
then compile
```python
cd build && make
```
To run do
```python
./pi0_evt_selection_study
```
with the output Root file in `./build/out.root`

---

To add histograms open `hists.json` and follow the examples for either 1D
```json
    {
      "name" : "<histo_name 1D>",
      "type": "TH1D",
      "axes": "<Title>>;<Xaxis>>;<Yaxis>>",
      "bins": 50,
      "u_lim": 200,
      "l_lim": -300
    }
```
or 2D histograms
```json
    {
      "name" : "<histo_name 2D>",
      "type": "TH2I",
      "axes": "<Title>>;<Xaxis>>;<Yaxis>>",
      "xbins": 10,
      "ybins": 10,
      "xu_lim": 10,
      "xl_lim": 0,
      "yu_lim": 10,
      "yl_lim": 0
    }
```
and add into the code `Histogram_object.th1_hists["name_of_histogram"]->Fill( variable )`


---

### Selection Parameters

1. Beam ToF cut
    * Upper 90 [ns]
    * Lower 80 [ns] 

To select beam pi+ particles there is an upper/lower cut placed on the 
the ToF.

2. Beam Type Cut 
    * beam type == 13
    
Pandora tags the beam as either 13 (track) or 11 (shower).
Since we only consider a pi+ primary we are looking for a track. 