BEGIN{echoline=1;}
/\/\*CUT HERE\*\// {echoline=0;}
/.*/               {if(echoline) print $0;}
