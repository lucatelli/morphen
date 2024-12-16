```
      %&&&&&+   +&&&&&&&*+                                                                 
      #@@@@@+  *#@@@@@#                                                                 *+ 
     *@@#*@@% *@@@%%@@+   %#&&@#%   &##&&#&%  +#@#&&##%+ +&&    &@@* &###&%%%*  +%&+   &@% 
    +@@@% &@@%@@#+ #@%  %@@*   @@& +@@%  +@@%  @@&  +@@% &@%   +@@%  &@#+      +@@@&  *@#  
    %@@#  +@@@@#  %@@  %@@*    @@& &@@&&&#&%  %@@%*%@#% &@@##&&@@@  +@@#&&%*+  #@#&@& #@*  
   +@@@*   &@@#   #@%  *@@%   &@#++@@%**#@%   @@&**+    &@*   %@@%  &@@+      &@@  *@#@&   
    %@&    *@&   &#&    +&##&#&%  &#&   *@@% %##+      *#%    ##&  %####&%%%*%@@&   *##+   
            +    +                       +%@&%***                          +&@@&+          
                     ++**%%%%%%%***++      +***+                         +%@@#%            
             +*%%&#@@@#&&%%%%%%%%&&#@@###@@@@@@#*                      *&@@#*              
         +*&#@@@#&%*+                *@@@@#&&@@@@                  +%#@@#%+                
      +&#@@@&*+                      #@@@*+  #@@@&%*++    +++**%&#@@#&*                    
    +&@@@&*                          &@@@@@@@@@@%%&####@@@@@##&&%*+                        
   %@@@%+                             +%&&&##&%+                                           
  %@@#+                                                                                    
  %&+                                                                                   
```
Older development version of `morphen` was moved to https://github.com/lucatelli/morphen_dev.

***Readme file under development!***

## What is `morphen`?
`morphen` puts together a collection of Python-based astronomical functionalities for image 
analysis and processing. 

The development version of `morphen` is available at: [morphen_dev](https://github.com/lucatelli/morphen_dev)

A major change in the code base is underway, and each one of the 
functionalities in the repository aforementioned will be migrated to individual repositories soon.
In the meantime, all tasks can be executed through this repository.


## Getting Started
### Installation
Currently, there is no option to install `morphen` (via `pip` or `conda`). 
However, installation and usage are straightforward. The code can be used as a module, 
interactively via Jupyter notebooks,
or via the command line interface (see "Important notes" below). For now, we recommend  
using it via Jupyter notebooks using the set of provided examples.

To install `morphen`, clone the repository and install using conda through the `environment.yml` file:
```bash
git clone https://github.com/lucatelli/morphen.git
cd morphen
conda env create -f environment.yml
conda activate morphen
```

See ([docs/install_instructions.md](docs/install_instructions.md)) for manual 
step-by-step installation instructions. There is also a experimental (manual) option
to install `morphen` on Mac OS with `arm64` architecture (M* chips).
An experimental `environment_mac.yml` file is available for this purpose. In principle (once 
`rosetta` and `miniconda` are installed), the creation of the `morphen` environment should be 
symply done by:
```commandline
arch -x86_64 conda env create -f environment_mac.yml
```






