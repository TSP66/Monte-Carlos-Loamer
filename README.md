# Monte Carlos Loamer

> __Loaming: "A method of geochemical prospecting in which samples of soil or other surficial material are tested for traces of the metal desired, its presence presumably indicating a near-surface orebody." - The American Geosciences Institute__

This is a Julia tool that employs [Monte Carlo methods](https://en.wikipedia.org/wiki/Monte_Carlo_method) to simulate the spread of geochemical anomalies mineral deposits could present on mountain sides. It is desgined with the gold prospector in mind - whereby these 'geochemical anomalies' take the form of gold-specks. It is able to simulate loaming (samples to likely deposit) and reverse-loaming (deposit to likely samples), maximising the theoretical chances of finding a 'patch' using what very limited information may be avaliable. It is currently able to model 'point' outcrops and [lode/vein](https://en.wikipedia.org/wiki/Lode) outcrops (including dip, strike and length). 

This tool could be of use in arid areas with a relatively thin overburden layer. The user is cautioned against using the tool in areas with significant overburden, glacial activity and/or landslides.

Functionality is being added to, suggestions are welcome.

An example of the log-anomaly presented by an auriferous reef on the slopes of the [Jukes-Darwin Mining Field](https://en.wikipedia.org/wiki/Mount_Jukes_mine_sites) from a variety of perspectives:

<table>
  <tr>
    <td>
      <img src="https://github.com/TSP66/Monte-Carlos-Loamer/blob/main/Images/Example_1.png" alt="Image 1">
    </td>
    <td>
      <img src="https://github.com/TSP66/Monte-Carlos-Loamer/blob/main/Images/Example_2.png" alt="Image 2">
    </td>
  </tr>
  <tr>
    <td>
      <img src="https://github.com/TSP66/Monte-Carlos-Loamer/blob/main/Images/Example_3.png" alt="Image 3">
    </td>
    <td>
      <img src="https://github.com/TSP66/Monte-Carlos-Loamer/blob/main/Images/Example_4.png" alt="Image 4">
    </td>
  </tr>
</table>

To make (useless and computationally expensive) colourful diagrams like this simply adjust the configurations in the config.toml file (instructions are in comments in the file) and run the sim.jl file. Elevation data is widely avaliable (in Australia) from [Elivs](https://elevation.fsdf.org.au/).

PlotlyJS, Images, ImageView, ArchGDAL, StatsBase and DelimitedFiles will need to be installed prior to running.

### Some basic performance notes:
- Keep either _res_ small (less than ~1000 if possible) or use smaller TIF files as PlotlyJS struggles to display larger images and simulation time also drastically increases with larger values of _res_
- It is possible to produce very clean images with not that many simulations (~1000 per sample), unless _display.log_ is true, in which case ~10,000 simulations per sample is needed. Simulating deposits with _display.log_ true can require in the order of 100,000s of simulations to produce a crisp images. Simulating deposits also tends to demand a higher _res_ value