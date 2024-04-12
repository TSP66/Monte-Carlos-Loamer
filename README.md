# Monte Carlos Loamer

> __Loaming: "A method of geochemical prospecting in which samples of soil or other surficial material are tested for traces of the metal desired, its presence presumably indicating a near-surface orebody." - The American Geosciences Institute__

This is a Julia tool that employs [Monte Carlo methods](https://en.wikipedia.org/wiki/Monte_Carlo_method) to simulate the spread of geochemical anomalies mineral deposits could present on mountain sides. It is desgined with the gold prospector in mind - whereby these 'geochemical anomalies' take the form of gold-specks. It is able to simulate loaming (samples to likely deposit) and reverse-loaming (deposit to likely samples), maximising the theoretical chances of finding a 'patch' using what very limited information may be avaliable. It is currently able to model 'point' outcrops and [lode/vein](https://en.wikipedia.org/wiki/Lode) outcrops (including dip, strike and length). 

This tool could be of use in arid areas with a relatively thin overburden layer. The user is cautioned against using the tool in areas with significant overburden, glacial activity and/or landslides.

Functionality is being added to, suggestions are welcome.

An example of an auriferous reef on the slopes of the [Jukes-Darwin Mining Field](https://en.wikipedia.org/wiki/Mount_Jukes_mine_sites):
![An example of a auriferous reef](https://github.com/TSP66/Monte-Carlos-Loamer/blob/main/Images/Example_1.png)

<table>
  <tr>
    <td>
      <img src="https://github.com/TSP66/Monte-Carlos-Loamer/blob/main/Images/Example_1.png" alt="Image 1">
    </td>
    <td>
      <img src="https://github.com/TSP66/Monte-Carlos-Loamer/blob/main/Images/Example_1.png" alt="Image 2">
    </td>
  </tr>
  <tr>
    <td>
      <img src="https://github.com/TSP66/Monte-Carlos-Loamer/blob/main/Images/Example_1.png" alt="Image 3">
    </td>
    <td>
      <img src="https://github.com/TSP66/Monte-Carlos-Loamer/blob/main/Images/Example_1.png" alt="Image 4">
    </td>
  </tr>
</table>

