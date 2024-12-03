theme_row_mid <- function(){
  theme_classic() %+replace%
    theme(axis.text = element_text(size = 5),
          axis.text.x = element_text(angle = 45, hjust = 0.5),
          axis.title = element_text(size = 5),
          axis.ticks.y = element_blank(),
          axis.title.y= element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank()) 
      
}

theme_row_border <- function(){
  theme_classic() %+replace%
    theme(axis.text = element_text(size = 5),
          axis.text.x = element_text(angle = 45, hjust = 0.5),
          axis.title = element_text(size = 5))
  
}

theme_col_mid <- function(){
  theme_classic() %+replace%
    theme(axis.text = element_text(size = 5),
          axis.title.y = element_text(size = 5),
          axis.text.y = element_text(angle = 45, vjust = 0.5),
          axis.ticks.x = element_blank(),
          axis.title.x= element_blank(),
          axis.text.x = element_blank(),
          axis.line.x = element_blank()) 
  
}

theme_col_border <- function(){
  theme_classic() %+replace%
    theme(axis.text = element_text(size = 5),
          axis.title.y = element_text(size = 5),
          axis.text.y = element_text(angle = 45, vjust = 0.5),
          axis.title = element_text(size = 5))
}

