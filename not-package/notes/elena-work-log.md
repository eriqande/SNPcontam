# Elena's Worklog
## 5/27/14
### Setting up Git Repository
* Set up a Github account under the username "evenable".
* Set up a folder "Hollings" in Terminal
  * `cd` moves to home directory
  * `mkdir Hollings` makes Hollings directory
* Cloned the SNPcontam repository using the url https://github.com/eriqande/SNPcontam.git.
  * `git clone https://github.com/eriqande/SNPcontam.git` clones the SNPcontam repository to my Hollings directory
* Changed some of the git configurations
  * `git config --global user.name "Elena Venable"` to input name
  * `git config --global user.email Elena_venable@brown.edu`: to input emai address
  * `git config --global core.editor emacs` changed the default editor to emacs
* Other important git commands
  * `git status` 
  * `git push origin master` pushes my local changes to SNPcontam repository
* Trying to get password caching to work
  * We tried this: `git config --global credential.helper osxkeychain`

### Work in Rstudio
* Downloaded packages needed for SNPcontam project
  * `knitr`, `devtools`, and `roxygen2`
* Git can also be accessed and used from Rstudio undert he "Git" tab
  

