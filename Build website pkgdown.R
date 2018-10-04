
 # From website: https://www.r-bloggers.com/building-a-website-with-pkgdown-a-short-guide/

 require(devtools)
 use_readme_rmd()
 use_news_md()
 #use_vignette("test")  #substitute with the name of your package
 use_vignette("qgg")
 
 use_github_links()
 
 #use_travis() : not sure we want to use this function. See devtools use_travis()
 
 #Create bare-bones website
 
 
 
 #Once your .yaml file is complete, just run
 build_site()
 #again and check the results. Then iterate ad-libitum until you are satisfied by the resulting structure.
 
 