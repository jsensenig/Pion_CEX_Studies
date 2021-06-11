###Code to study vertex reconstruction.

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
./reco_vtx_study
```
with the output Root file `./build/out.root`

---

To add histograms open `hists.json` and follow the examples for either 1D
```json
    {
      "name" : "hEndXDiffEl",
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
      "name" : "Example 2D Histo",
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
and add the line into the code as `Histogram_object.th1_hists["name_of_histogram"]->Fill( variable )`

---

Event selection criteria

1. Beam type (track = 13 or shower = 11)
2. ToF cut
3. IsBeamLike (BI to TPC position match + angle)
4. 2 daughter tracks (pi+ and >0 proton) no showers using track/shower CNN score
5. 1 pi+ and >0 proton daughters using dE/dx vs Residual range fitting PID