import numpy as np  
import matplotlib.pyplot as plt
import os 
from string import Template
import subprocess

def setIPG( ipg ):   
   ipg = int(ipg)
   if (ipg == 1) : 
      method = 'NIPG'
   elif (ipg == -1): 
      method = 'SIPG'
   elif (ipg == 0): 
      method = 'IIPG' 
   else:
      method = 'Unknown method'
   return method

# load the errs and other parameters of the computation from the  file 
def fromFileToLists(directory): 
   with open(directory) as f: 
      case, eps, ipg = \
      [float(x) for x in next(f).split()] # read the first line 
      line = [float(x) for x in next(f).split()]
      line2 = [float(x) for x in next(f).split()] 
      f.close() 
   
   return line, line2, case, eps, ipg 

# make tex Table from the given files 
def makeTex(directories):  
      
      lines = [] 
      lines2 = [] 
      cases = [] 
      epsilons = []
      ipgs = []
      for directory in directories: 
         line, line2, case, eps, ipg = fromFileToLists(directory) 
         lines.append(line)
         lines2.append(line2) 
         cases.append(case)
         epsilons.append(eps) 
         ipgs.append(ipg) 
      print 'Lines type:', type(lines)
      print 'Lines:', lines   
      
      # CONTROLS- we work the only with case, eps, ipg
      for c in cases: 
         if (c != cases[0]): 
            print 'The case are not the same!'
            quit()             
      for c in epsilons: 
         if (c != epsilons[0]): 
            print 'Epsilons are not the same!'
            quit() 
      for c in ipgs: 
         if (c != ipgs[0]): 
            print 'IPG are not the same!'
            quit() 
               
      
      file_name = 'output_case_' + str(case)
      file_name_tex = file_name + '.tex'
      #if not os.path.isdir( paramets.path ): os.makedirs(paramets.path)
      #full_path = os.path.join(paramets.path, file_name)  
      
      n_lines = len( lines ) 
      
      mIpg = setIPG(ipg)
      
      beginning =  Template(r'''
\documentclass[10pt]{article}
\usepackage{epsfig}
\usepackage{rotating}
\usepackage{amssymb}
\usepackage{amsmath}        

% soucast baliku texlive-sciene
% nerozumi si s nekterymi definicemi napr \mA, je treba je zmenit ( mA se pouziva jako miliamper), atd.
\usepackage[scientific-notation=true]{siunitx}
\sisetup{round-mode=places,round-precision=1,tight-spacing = true}
\newcommand{\numeff}[1]{\num[round-mode=places,round-precision=2]{#1}}
\newcommand{\numm}[1]{\num[round-mode=places,round-precision=3]{#1}}
% \numeff - effectivity index format, \num 1.23 \times 10^{-3} format

\renewcommand{\arraystretch}{1.2}

\hoffset=-41mm
\voffset=-8mm
\topmargin -2cm
%\leftmargin -25mm

\begin{document} 
\textbf{\large RTNst estimates \\ } 
\\
\textbf{ Type of problem: $case \\} 
\textbf{ Epsilon = \numeff{$epsilon} \\ } 
\textbf{ type of DG = $ipg \\ } 

            ''')
      table_beg = r'''
{\scriptsize             %  mensi
%{\tiny                   %  nejmensi
\begin{tabular}{|ccccc|rrrr|rrr|}
\hline
$N_{h}$ & $h$ & $\tau$ & $p$ & $q$ &  
 $||\psi||_{Z}$ & $\eta$ & $ ||e_{FR}||$ & $\eta_{NC}$ &
$i_{eff}$ & $i_{eff}^{tot}$ & $i_{eff}^{fr}$ \\ 
\hline           
         ''' 
      table2_beg = r'''
{\scriptsize             %  mensi
%{\tiny                   %  nejmensi
\begin{tabular}{ccccc|rrrr}
\hline
$N_{h}$ & $h$ & $\tau$ & $p$ & $q$ &  
$\eta_R$ & $\eta_F$ & $\eta_T$ & $\eta_{NC} $ \\ 
\hline           
         ''' 
         
      table_end = r''' 
\hline
\end{tabular}
} %%% end of font size  
%\end{document}
      '''
      
      doc_end = r''' 
\end{document}
      '''

      first = beginning.substitute( case = case, \
                          epsilon = eps, ipg = mIpg )
         
      with open( file_name, 'w') as f:
         f.write( first )    
         # First table 
         f.write( table_beg )  
         for i in range( n_lines ): 
            f.write( '{:5} {:2} {:8}{:12.12f}{:1} {:2} {:8}{:12.12f}{:1} {:2} {:5.0f} {:2} {:5.0f}'.format( \
                     int(lines[i][0]) , '&', '\\numeff{', lines[i][1], '}' , '&', '\\numeff{',  lines[i][2], '}', '&', lines[i][3] , '&', lines[i][4] ) )           
            # errors
            for j in range(5,9): 
               f.write( '{:2} {:6} {:12.12f}{:1}'. format( '&' , '\\numm{' , lines[i][j] , '}' ) ) 
            # i_eff
            for j in range(9,12): 
               f.write( '{:2}{:8}{:12.12f}{:1} '. format( '&', '\\numeff{', lines[i][j] , '}' ) ) 
            f.write(  '{:2}'.format('\\\\ \n') )
     
         f.write( table_end )     
         # SECOND table 
         f.write( table2_beg )  
         for i in range( n_lines ): 
            f.write( '{:5} {:2} {:8}{:12.12f}{:1} {:2} {:8}{:12.12f}{:1} {:2} {:5.0f} {:2} {:5.0f}'.format( \
                     int(lines[i][0]) , '&', '\\numeff{', lines[i][1], '}' , '&', '\\numeff{',  lines[i][2], '}', '&', lines[i][3] , '&', lines[i][4] ) )           
            for j in range(4): 
               f.write( '{:2} {:6} {:12.12f}{:1}'. format( '&' , '\\numm{' , lines2[i][j] , '}' ) )
            f.write(  '{:2}'.format('\\\\ \n') ) 
         
         f.write( table_end )     
         f.write( doc_end)
         f.close()         
      
      cmd = ['pdflatex', '-interaction', 'nonstopmode', file_name]
      proc = subprocess.Popen(cmd)
      proc.communicate()

      retcode = proc.returncode
      if not retcode == 0:
         os.unlink( file_name + '.pdf')
         raise ValueError('Error {} executing command: {}'.format(retcode, ' '.join(cmd))) 
      #full_pdf = os.path.join(paramets.path, "output.pdf")  
      #os.rename( "output.pdf" , full_pdf )
      
      # I do not want to remove tex file
      # os.unlink(  os.path.join(paramets.path, 'output.tex') )
      os.unlink( file_name + '.log') 
      os.unlink( file_name + '.aux')
      
      print 'The file', file_name, '.pdf was created!'




