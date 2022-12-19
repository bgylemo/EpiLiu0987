# Made by to Antonio Lentini (https://github.com/ALentini)

require(ggplot2)
theme_AL_simple <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(colour="black"),axis.line=element_line(colour="black"), ... )
}
theme_AL_simple_rotX <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line=element_line(colour="black"),axis.text = element_text(colour="black"),axis.text.x=element_text(angle = 90,hjust=1), ... )
}
theme_AL_box <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(colour="black"),axis.line=element_line(colour="black"), ... )
}
theme_AL_box_rotX <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(colour="black"),axis.line=element_line(colour="black"), axis.text.x=element_text(angle = 90,hjust=1), ... )
}
theme_AL_box_rotX_45 <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(colour="black"),axis.line=element_line(colour="black"), axis.text.x=element_text(angle = 45), ... )
}
theme_AL_simple_rotX_45 <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line=element_line(colour="black"),axis.text = element_text(colour="black"),axis.text.x=element_text(angle = 45), ... )
}