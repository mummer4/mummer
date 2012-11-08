package TIGR::Foundation;
{
 
=head1 NAME

TIGR::Foundation - TIGR Foundation object

=head1 SYNOPSIS

  use TIGR::Foundation;
  my $obj_instance = new TIGR::Foundation;

=head1 DESCRIPTION

This module defines a structure for Perl programs to utilize
logging, version reporting, and dependency checking in a simple way.

=cut

   BEGIN {
      require 5.006_00;                       # error if using Perl < v5.6.0  
   }

   use strict;
   use Cwd;
   use Cwd 'chdir';
   use Cwd 'abs_path';
   use File::Basename;
   use Getopt::Long;
   use IO::Handle;
   use POSIX qw(strftime);
   use Sys::Hostname;
   use English;   

   require Exporter;

   our @ISA;
   our @EXPORT;
   @ISA = ('Exporter');
   @EXPORT = qw(
                isReadableFile
                isWritableFile
                isExecutableFile
                isCreatableFile
                isReadableDir
                isWritableDir
                isCreatableDir
                isCreatablePath

                getISODate
                getSybaseDate
                getMySQLDate
                getFilelabelDate
                getLogfileDate
               );

   ## internal variables and identifiers
   our $REVISION = (qw$Revision: 1.1.1.1 $)[-1];
   our $VERSION = '1.1'; 
   our $VERSION_STRING = "$VERSION (Build $REVISION)";
   our @DEPEND = ();                          # there are no dependencies


   ## prototypes

   # Functional Class : general
   sub new();
   sub getProgramInfo($);
   sub runCommand($);

   # Functional Class : depend
   sub printDependInfo();
   sub printDependInfoAndExit();
   sub addDependInfo(@);

   # Functional Class : version
   sub getVersionInfo();
   sub printVersionInfo();
   sub printVersionInfoAndExit();
   sub setVersionInfo($);

   # Functional Class : help
   sub printHelpInfo();
   sub printHelpInfoAndExit();
   sub setHelpInfo($);

   # Functional Class : usage
   sub printUsageInfo();
   sub printUsageInfoAndExit();
   sub setUsageInfo($);

   # Functional Class : files
   sub isReadableFile($);
   sub isExecutableFile($);
   sub isWritableFile($);
   sub isCreatableFile($);
   sub isReadableDir($);
   sub isWritableDir($);
   sub isCreatableDir($);
   sub isCreatablePath($);

   # Functional Class : date
   sub getISODate(;@);
   sub getSybaseDate(;@);
   sub getMySQLDate(;@);
   sub getFilelabelDate(;@);
   sub getLogfileDate(;@);

   # Functional Class : logging
   sub setDebugLevel($;$);
   sub getDebugLevel();
   sub setLogFile($;$);
   sub getLogFile();
   sub getErrorFile();
   sub printDependInfo();
   sub invalidateLogFILES();
   sub cleanLogFILES();
   sub closeLogERROR();
   sub closeLogMSG();
   sub openLogERROR();
   sub openLogMSG();
   sub logAppend($;$);
   sub debugPush();
   sub debugPop();
   sub logLocal($$);
   sub logError($;$);
   sub bail($;$);

   # Functional Class : modified methods
   sub TIGR_GetOptions(@);

   ## Implementation
   

# Functional Class : general

=over

=item $obj_instance = new TIGR::Foundation;

This function creates a new instance of the TIGR::Foundation
object.  A reference pointing to the object is returned on success.  Otherwise,
this method returns undefined.

=cut


   sub new() {
      
      my $self = {};
      my $pkg = shift;
      my $user_name = getpwuid($<);
      my $host_name = hostname();

      # create the object
      bless $self, $pkg;

      ## Instance variables and identifiers, by functional class
      
      # Functional Class : general
      my $pname  = basename($0, ()); # extract the script name
       
      if((defined ($pname)) && ($pname =~ /^(.*)$/)) {
          $pname = $1;
	  $self->{program_name} = $pname ;
      }
      
      if ($self->{program_name} =~ /^-$/) {     # check if '-' is the input
         $self->{program_name} = "STDIN";
      }
      
      my $pcommand = join (' ', @ARGV);
      
      if((defined ($pcommand)) && ($pcommand =~ /^(.*)$/)) {
          $pcommand = $1;
	  $self->{invocation} = $pcommand ;
      }
      
      # The following four variables are to contain information specified by
      # the 'host' program; there are methods of setting and retrieving each

      # Functional Class : depend
      @{$self->{depend_info}} = ();
   
      # Functional Class : version
      $self->{version_info} = undef; 
   
      # Functional Class : help
      $self->{help_info} = undef; 
   
      # Functional Class : usage
      $self->{usage_info} = undef; 
   
      # Functional Class : logging
      $self->{debug_level} = undef;             # default debug is not defined
      @{$self->{debug_store}} = ();             # the backup debug level stack
      @{$self->{debug_queue}} = ();             # queue used by MSG routine
      @{$self->{error_queue}} = ();             # queue used by ERROR routine
      $self->{max_debug_queue_size} = 100;      # maximum size for queue before
                                                # log entries are expired
      @{$self->{log_files}} =                   # these log files are consulted
         ("$self->{program_name}.log",          # on file write error and are
          "/tmp/$self->{program_name}.$$.log"); # modified by setLogFile
      $self->{msg_file_open_flag} = 0;          # flag to check logLocal file
      $self->{error_file_open_flag} = 0;        # flag to check logError file
      $self->{msg_file_used} = 0;               # flag to indicate if log file
      $self->{error_file_used} = 0;             #   has been written to
      $self->{msg_append_flag} = 0;             # by default logs are truncated
      $self->{error_append_flag} = 0;           # by default logs are truncated
      $self->{log_append_setting} = 0;          # (truncate == 0)
      $self->{static_log_file} = undef;         # user defined log file
      $self->{start_time} = undef;              # program start time
      $self->{finish_time} = undef;             # program stop time

      # Log program invocation
      $self->logLocal("START: " . "Username:$user_name, ".
                      "Hostname: $host_name ". $self->getProgramInfo('name') .
                      " " . $self->getProgramInfo('invocation'), 0);
      $self->{start_time} = time;

      return $self;
   }



=item $value = $obj_instance->getProgramInfo($field_type);

This function returns field values for specified field types describing
attributes of the program.  The C<$field_type> parameter must be a listed
attribute: C<name>, C<invocation>, C<env_path>, C<abs_path>.
The C<name> field specifies the bare name of the executable.  The
C<invocation> field specifies the command line arguments passed to the
executable.   The C<env_path> value returns the environment path to the
working directory.  The C<abs_path> value specifies the absolute path to the
working directory.  If C<env_path> is found to be inconsistent, then that
value will return the C<abs_path> value.  If an invalid C<$field_type> is 
passed, the function returns undefined.  

=cut


   sub getProgramInfo($) {
      my $self = shift;
      my $field_type = shift;
      my $return_value = undef;
      if (defined $field_type) {
         $field_type =~ /^name$/ && do {
            $return_value = $self->{program_name};
         };
         $field_type =~ /^invocation$/ && do {
            $return_value = $self->{invocation};
         };
         $field_type =~ /^env_path$/ && do {
            my $return_value = "";
            if (
                (defined $ENV{'PWD'}) &&
                (abs_path($ENV{'PWD'}) eq abs_path(".")) &&
                ($ENV{'PWD'} =~ /^(.*)$/)     
               ) {
	       $ENV{'PWD'} = $1;
               $return_value = $ENV{'PWD'};
            }
            else {
	       my $tmp_val = abs_path(".");

               if((defined ($tmp_val)) && ($tmp_val =~ /^(.*)$/)) {
                  $tmp_val = $1;
	          $return_value = $tmp_val;
               }
            }
            return $return_value;
         };

         $field_type =~ /^abs_path$/ && do {
            my $tmp_val = abs_path(".");

            if((defined ($tmp_val)) && ($tmp_val =~ /^(.*)$/)) {
               $tmp_val = $1;
	       $return_value = $tmp_val;
            }
         };
      }
      return $return_value;
   }

=item $exit_code = $obj_instance->runCommand($command_str);

This function passes the argument C<$command_str> to /bin/sh
for processing.  The return value is the exit code of the 
C<$command_str>.  If the exit code is not defined, then either the signal or
core dump value of the execution is returned, whichever is applicable.  Perl
variables C<$?> and C<$!> are set accordingly.  If C<$command_str> is not 
defined, this function returns undefined.  Log messages are recorded at log
level 4 to indicate the type of exit status and the corresponding code.
Invalid commands return -1.

=cut


   sub runCommand($) {
       my $self = shift;
       my $command_str = shift;
       my $exit_code = undef; 
       my $signal_num = undef;
       my $dumped_core = undef;
       my $invalid_command = undef;
       my $return_value = undef;
       my @info_arr = getpwuid($<);
       my $len = @info_arr;
       my $home_dir = $info_arr[7];
       my $current_dir = $self->getProgramInfo("abs_path"); 
       
       if((defined ($ENV{PATH})) && ($ENV{PATH} =~ /^(.*)$/)) {#taint checking
	  $ENV{PATH} = $1; 
          my $path_var = $ENV{PATH};
          my @paths = split /:/, $path_var;
          my $pathval = undef;
          my $i = 0;
          my $paths_len = @paths;

          for ($i = 0; $i < $paths_len ; $i++) {
             #substituting ~ with the home pathname.
	     $pathval = $paths[$i];
	     $pathval =~ s/^~$/$home_dir/g;
             my $home_root = $home_dir."\/";
             $pathval =~ s/^~\//$home_root/g;
             
             #substituting . with the current pathname.
             $pathval =~ s/^\.$/$current_dir/g;
             my $current_root = $current_dir."\/";
             $pathval =~ s/^\.\//$current_root/g;
             $paths[$i] = $pathval;
	  }
          
          $ENV{PATH} = join(":", @paths);
       }

       if((defined ($command_str)) && ($command_str =~ /^(.*)$/)) {#taint 
                                                                   #checking
          $command_str = $1;
          system($command_str);
          $exit_code = $? >> 8;
          $signal_num = $? & 127;
          $dumped_core = $? & 128;

          if ($? == -1) {
             $invalid_command = -1;
          } 
 
          if ( 
             (!defined $invalid_command) &&
             ($exit_code == 0) &&
             ($signal_num == 0) &&
             ($dumped_core != 0)
             ) {
            
             $self->logLocal("Command '" . $command_str . "' core dumped", 4);
             $return_value = $dumped_core;
           }
          elsif (
                (!defined $invalid_command) &&
                ($exit_code == 0) &&
                ($signal_num != 0)
                ) {
             
             $self->logLocal("Command '" . $command_str .
                             "' exited on signal " . $signal_num, 4);
             $return_value = $signal_num;
          }
          elsif ((!defined $invalid_command)) {
             
             $self->logLocal("Command '" . $command_str .
                             "' exited with exit code " . $exit_code, 4);
             $return_value = $exit_code;
          }
          else {
           
             $self->logLocal("Command '" . $command_str .
                             "' exited with invalid code " . $?, 4);
             $return_value = $?;
          }
       }
       return $return_value;
   }
   

# Functional Class : depend

=item $obj_instance->printDependInfo();

The C<printDependInfo()> function prints the dependency list created by
C<addDependInfo()>.  One item is printed per line.

=cut


   sub printDependInfo() {
      my $self = shift;
      foreach my $dependent (@{$self->{depend_info}}) {
         print STDERR $dependent, "\n";
      }
   }


=item $obj_instance->printDependInfoAndExit();

The C<printDependInfoAndExit()> function prints the dependency list created by
C<addDependInfo()>.  One item is printed per line.  The function exits with
exit code 0. 

=cut


   sub printDependInfoAndExit() {
      my $self = shift;
      $self->printDependInfo();
      exit 0;
   }


=item $obj_instance->addDependInfo(@depend_list);

The C<addDependInfo()> function adds C<@depend_list> information
to the dependency list.  If C<@depend_list> is empty, the internal
dependency list is emptied.  Contents of C<@depend_list> are not checked
for validity (eg. they can be composed entirely of white space or
multiple files per record).  The first undefined record in C<@depend_list>
halts reading in of dependency information.

=cut


   sub addDependInfo(@) {
      my $self = shift;
      my $num_elts = 0;
      while (my $data_elt = shift @_) {
         push (@{$self->{depend_info}}, $data_elt);
         $num_elts++;
      }
      if ($num_elts == 0) {
         @{$self->{depend_info}} = ();
      }
   }


# Functional Class : version

=item $version_string = $obj_instance->getVersionInfo();

The C<getVersionInfo()> function returns the version information set by the
C<setVersionInfo()> function.

=cut


   sub getVersionInfo() {
      my $self = shift;
      return $self->{version_info};
   }


=item $obj_instance->printVersionInfo();

The C<printVersionInfo()> function prints the version information set by the
C<setVersionInfo()> function.  If there is no defined version information,
a message is returned notifying the user.

=cut


   sub printVersionInfo() {
      my $self = shift;
      if (defined $self->getVersionInfo()) {
         print STDERR $self->getProgramInfo('name'), 
            " ", $self->getVersionInfo(), "\n";
      }
      else {
         print STDERR $self->getProgramInfo('name'),
            " has no defined version information\n";
      }
   }


=item $obj_instance->printVersionInfoAndExit();

The C<printVersionInfoAndExit()> function prints version info set by the
C<setVersionInfo()> function.  If there is no defined version information,
a message is printed notifying the user.  This function calls exit with
exit code 0. 

=cut


   sub printVersionInfoAndExit() {
      my $self = shift;
      $self->printVersionInfo();
      exit 0;
   }


=item $obj_instance->setVersionInfo($version_string);

The C<setVersionInfo()> function sets the version information to be reported
by C<getVersionInfo()>.  If C<$version_string> is empty, invalid, or
undefined, the stored version information will be undefined.

=cut


   sub setVersionInfo($) {
      my $self = shift;
      my $v_info = shift;
      if (
          (defined $v_info) &&
          ($v_info =~ /\S/) &&
          ((ref $v_info) eq "")
         ) {
         $self->{version_info} = $v_info;
      }
      else {
         $self->{version_info} = undef;
      }
   }


# Functional Class : help

=item $obj_instance->printHelpInfo();

The C<printHelpInfo()> function prints the help information passed by the
C<setHelpInfo()> function.

=cut


   sub printHelpInfo() {
      my $self = shift;
      if (defined $self->{help_info}) {
         print STDERR $self->{help_info};
      }
      else {
         print STDERR "No help information defined.\n";
      }
   }


=item $obj_instance->printHelpInfoAndExit();

The C<printHelpInfoAndExit()> function prints the help info passed by the
C<setHelpInfo()> function.  This function exits with exit code 0.

=cut


   sub printHelpInfoAndExit() {
      my $self = shift;
      $self->printHelpInfo();
      exit 0;
   }


=item $obj_instance->setHelpInfo($help_string);

The C<setHelpInfo()> function sets the help information via C<$help_string>.
If C<$help_string> is undefined, invalid, or empty, the help information 
is undefined.

=cut


   sub setHelpInfo($) {
      my $self = shift;
      my $help_string = shift;
      if (
          (defined $help_string) &&
          ($help_string =~ /\S/) &&
          ((ref $help_string) eq "")
         ) {
	 chomp($help_string);#removing a new line if it is there.
         $self->{help_info} = $help_string."\n";#adding a new line to help.
      }
      else {
         $self->{help_info} = undef;
      }
   }


# Functional Class : usage

=item $obj_instance->printUsageInfo();

The C<printUsageInfo()> function prints the usage information reported by the
C<setUsageInfo()> function.  If no usage information is defined, but help
information is defined, help information will be printed.

=cut


   sub printUsageInfo() {
     
      my $self = shift;
      if (defined $self->{usage_info}) {
         print STDERR $self->{usage_info};
      }
      elsif (defined $self->{help_info}) {
         print STDERR $self->{help_info};
      }
      else {
         print STDERR "No usage information defined.\n";
      }
   }


=item $obj_instance->printUsageInfoAndExit();

The C<printUsageInfoAndExit()> function prints the usage information the
reported by the C<setUsageInfo()> function and exits with status 1.

=cut


   sub printUsageInfoAndExit() {
      my $self = shift;
      $self->printUsageInfo();
      $self->bail("Incorrect command line");
   }


=item $obj_instance->setUsageInfo($usage_string);

The C<setUsageInfo()> function sets the usage information via C<$usage_string>.
If C<$usage_string> is undefined, invalid, or empty, the usage information 
is undefined.

=cut


   sub setUsageInfo($) {
      my $self = shift;
      my $usage_string = shift;
      if (
          (defined $usage_string) &&
          ($usage_string =~ /\S/) &&
          ((ref $usage_string) eq "")
         ) {
	 chomp($usage_string); #removing a new line if it is there.
         $self->{usage_info} = $usage_string."\n";#adding a new line to usage
      }
      else {
         $self->{usage_info} = undef;
      }
   }


# Functional Class : files

=item $valid = isReadableFile($file_name);

This function accepts a single scalar parameter containing a file name.
If the file corresponding to the file name is a readable plain file or symbolic
link, this function returns 1.  Otherwise, the function returns 0.  If the file
name passed is undefined, this function returns 0 as well.

=cut


   sub isReadableFile($) {
      my $file = shift;
      if (scalar(@_) != 0) { #incase the method was invoked as an instance 
                             #method
         $file  = shift;
      }
 
      if (defined ($file) &&             # was a file name passed?
          ((-f $file) || (-l $file)) &&  # is the file a file or sym. link?
          (-r $file)                     # is the file readable?
         ) {
         return 1;
      }
      else {
         return 0;
      }
   }


=item $valid = isExecutableFile($file_name);

This function accepts a single scalar parameter containing a file name.
If the file corresponding to the file name is an executable plain file
or symbolic link, this function returns 1.  Otherwise, the function returns 0.
If the file name passed is undefined, this function returns 0 as well.

=cut


   sub isExecutableFile($) {
      my $file = shift;
       if (scalar(@_) != 0) { # incase the method was invoked as a instance
                              # method
         $file  = shift;
      }
 
      if (defined ($file) &&             # was a file name passed?
          ((-f $file) || (-l $file)) &&  # is the file a file or sym. link?
          (-x $file)                     # is the file executable?
         ) {
         return 1;
      }
      else {
         return 0;
      }
   }


=item $valid = isWritableFile($file_name);

This function accepts a single scalar parameter containing a file name.
If the file corresponding to the file name is a writable plain file
or symbolic link, this function returns 1.  Otherwise, the function returns 0.
If the file name passed is undefined, this function returns 0 as well.

=cut


   sub isWritableFile($) {
      my $file = shift;
      if (scalar(@_) != 0) { # incase the method was invoked as a instance
                             #  method
         $file  = shift;
      }
 
      if (defined ($file) &&             # was a file name passed?
          ((-f $file) || (-l $file)) &&  # is the file a file or sym. link?
          (-w $file)                     # is the file writable?
         ) {
         return 1;
      }
      else {
         return 0;
      }
   }


=item $valid = isCreatableFile($file_name);

This function accepts a single scalar parameter containing a file name.  If
the file corresponding to the file name is creatable this function returns 1.
The function checks if the location of the file is writable by the effective
user id (EUID).  If the file location does not exist or the location is not
writable, the function returns 0.  If the file name passed is undefined,
this function returns 0 as well.  Note that files with suffix F</> are not
supported under UNIX platforms, and will return 0.

=cut


   sub isCreatableFile($) {
      my $file = shift;
      if (scalar(@_) != 0) { # incase the method was invoked as an instance
                             # method
         $file  = shift;
      }

      my $return_code = 0;

      if (
          (defined ($file)) &&
          (! -e $file) &&
          ($file !~ /\/$/) 
         ) {
         my $dirname = dirname($file);
         # check the writability of the directory
         $return_code = isWritableDir($dirname);
      }
      else {
         # the file exists, it's not creatable
         $return_code = 0;
      }
      return $return_code;
   }


=item $valid = isReadableDir($directory_name);

This function accepts a single scalar parameter containing a directory name.
If the name corresponding to the directory is a readable, searchable directory 
entry, this function returns 1.  Otherwise, the function returns 0.  If the
name passed is undefined, this function returns 0 as well.

=cut


   sub isReadableDir($) {
      my $file = shift;
      if (scalar(@_) != 0) { # incase the method was invoked as an instance
                             # method
         $file  = shift;
      }
 
      if (defined ($file) &&             # was a name passed?
          (-d $file) &&                  # is the name a directory?
          (-r $file) &&                  # is the directory readable?
          (-x $file)                     # is the directory searchable?
         ) {
         return 1;
      }
      else {
         return 0;
      }
   }


=item $valid = isWritableDir($directory_name);

This function accepts a single scalar parameter containing a directory name.
If the name corresponding to the directory is a writable, searchable directory 
entry, this function returns 1.  Otherwise, the function returns 0.  If the
name passed is undefined, this function returns 0 as well.

=cut


   sub isWritableDir($) {
      my $file = shift;
      if (scalar(@_) != 0) { # incase the method was invoked as an instance
                             # method
         $file  = shift;
      }
 
      if (defined ($file) &&             # was a name passed?
          (-d $file) &&                  # is the name a directory?
          (-w $file) &&                  # is the directory writable?
          (-x $file)                     # is the directory searchable?
         ) {
         return 1;
      }
      else {
         return 0;
      }
   }


=item $valid = isCreatableDir($directory_name);

This function accepts a single scalar parameter containing a directory name.  
If the name corresponding to the directory is creatable this function returns 
1. The function checks if the immediate parent of the directory is writable by
the effective user id (EUID).  If the parent directory does not exist or the 
tree is not writable, the function returns 0.  If the directory name passed is
undefined, this function returns 0 as well.

=cut


   sub isCreatableDir($) {
      my $dir = shift;
      if (scalar(@_) != 0) { # incase the method was invoked as an instance
                             # method
         $dir  = shift;
      }
      my $return_code = 0;

      if (defined ($dir)) {
         $dir =~ s/\/$//g;
         $return_code = isCreatableFile($dir);
      }
      return $return_code;
   }


=item $valid = isCreatablePath($path_name);

This function accepts a single scalar parameter containing a path name.  If
the C<$path_name> is creatable this function returns 1. The function checks 
if the directory hierarchy of the path is creatable or writable by the
effective user id (EUID).  This function calls itself recursively until
an existing directory node is found.  If that node is writable, ie. the path
can be created in it, then this function returns 1.  Otherwise, the function
returns 0.  This function also returns zero if the C<$path_name> supplied
is disconnected from a reachable directory tree on the file system.
If the path already exists, this function returns 0.  The C<$path_name> may
imply either a path to a file or a directory.  Path names may be relative or
absolute paths.  Any unresolvable relative paths will return 0 as well.  This
includes paths with F<..> back references to nonexistent directories.
This function is recursive whereas C<isCreatableFile()> and 
C<isCreatableDir()> are not.

=cut


   sub isCreatablePath($) {
      my $pathname = shift;
      if (scalar(@_) != 0) { # incase the method was invoked as an instance
                             # method
         $pathname  = shift;
      }
      my $return_code = 0;

      if (defined $pathname) {
         # strip trailing '/'
         $pathname =~ s/(.+)\/$/$1/g;
         my $filename = basename($pathname);
         my $dirname = dirname($pathname);
         if (
             (! -e $pathname) &&
             ($dirname ne $pathname) &&
             ($filename ne "..")
            ) {
            if (-e $dirname) {
               $return_code = isWritableDir($dirname);
            }
            else {
               $return_code = isCreatablePath($dirname);
            }
         }
         else {
            $return_code = 0;
         }
      }
      return $return_code;
   }
         

# Functional Class : date

=item $date_string = getISODate($tm);

This function returns the ISO 8601 datetime as a string given a time
structure as returned by the C<time> function.  If no arguments
are supplied, this function returns the current time.  If incorrect
arguments are supplied, this function returns undefined.

=cut


   sub getISODate(;@) {
      #checking if the function is invoked as an instance method.
      if((defined(ref $_[0])) && ((ref $_[0]) eq "TIGR::Foundation")){
         shift;
      } 
      my @time_val = @_;
      my $time_str = undef;
      if (scalar(@time_val) == 0) {
         @time_val = localtime;
      }
      eval {
         $time_str = strftime "%Y-%m-%d %H:%M:%S", @time_val;
      };
      return $time_str;
   }


=item $date_string = getSybaseDate(@tm);

This function returns a Sybase formatted datetime as a string given a time
structure as returned by the C<time> function.  If no arguments
are supplied, this function returns the current time.  If incorrect
arguments are supplied, this function returns undefined.  The date string
returned is quoted according to Sybase requirements.

=cut


   sub getSybaseDate(;@) {
      #checking if the function is invoked as an instance method.
      if((defined(ref $_[0])) && ((ref $_[0]) eq "TIGR::Foundation")){
         shift;
      } 
      my @time_val = @_;
      my $time_str = undef;
      if (scalar(@time_val) == 0) {
         @time_val = localtime;
      }
      eval {
         $time_str = strftime "\'%b %d %Y %I:%M%p\'", @time_val;
      };
      return $time_str;
   }
      

=item $date_string = getMySQLDate(@tm);

This function returns a MySQL formatted datetime as a string given a time
structure as returned by the C<time> function.  If no arguments
are supplied, this function returns the current time.  If incorrect
arguments are supplied, this function returns undefined.  The datetime string
returned is prequoted according to MySQL requirements.

=cut


   sub getMySQLDate(;@) {
      #checking if the function is invoked as an instance method.
      if((defined(ref $_[0])) && ((ref $_[0]) eq "TIGR::Foundation")){
         shift;
      }
      my @time_val = @_;
      my $time_str = undef;
      if (scalar(@time_val) == 0) {
         @time_val = localtime;
      }
      $time_str = getISODate(@time_val);
      if (defined $time_str) {
         $time_str = "\'$time_str\'";
      }
      return $time_str;
   }
      

=item $date_string = getFilelabelDate(@tm);

This function returns the date (not time) as a compressed string
suitable for use as part of a file name.  The format is YYMMDD.
The optional parameter should be a time structure as returned by 
the C<time> function.  If no arguments are supplied, the current time
is used.  If incorrect arguments are supplied, this function returns
undefined.

=cut


   sub getFilelabelDate(;@) {
      #checking if the function is invoked as an instance method.
      if((defined(ref $_[0])) && ((ref $_[0]) eq "TIGR::Foundation")){
         shift;
      }
      my @time_val = @_;
      my $time_str = undef;
      if (scalar(@time_val) == 0) {
         @time_val = localtime;
      }
      eval {
         $time_str = strftime "%y%m%d", @time_val;
      };
      return $time_str;
   }
      

=item $date_string = $obj_instance->getLogfileDate(@tm);

This function returns the datetime as a formatted string
suitable for use as a log entry header.  The optional parameter
should be a time structure as returned by the C<time> function.
If no arguments are supplied, this function uses the current time.
If incorrect arguments are supplied, this function sets the date/time fields
of the log entry string to C< INVALID|XXXXXX|>.

=cut


   sub getLogfileDate(;@) {
      #checking if the function is invoked as an instance method.
      if((defined(ref $_[0])) && ((ref $_[0]) eq "TIGR::Foundation")){
         shift;
      } 
      my @time_val = @_;
      my $time_str = undef;
      my $log_form = undef;
      if (scalar(@time_val) == 0) {
         @time_val = localtime;
      }
      eval {
         $time_str = strftime("%Y%m%d|%H%M%S|", @time_val);
      };
      if (!defined $time_str) {
         $time_str = " INVALID|XXXXXX|";
      }
      $log_form = $time_str . sprintf("%6d| ", $$);
      return $log_form;
   }
      

# Functional Class : logging

=item $obj_instance->setDebugLevel($new_level);

This function sets the level of debug reporting according to C<$new_level>.
If the debug level is less than 0, all debug reporting is turned off.
It is impossible to turn off error reporting from C<bail()>.  If C<$new_level>
is undefined, the debug level is set to 0.  This function maintains 
compatibility with C<GetOptions()>, and will accept a second parameter
the debug level, provided it is an integer.  In such cases, the first parameter
is checked only if the second parameter is invalid.  By default, the default
level is undefined.  To turn on debugging, you must invoke this function.

=cut


   sub setDebugLevel($;$) {
      my $self = shift;
      my $new_level = shift;
      my $getopts_new_level = shift;

      if (
          (defined $getopts_new_level) &&
          ($getopts_new_level =~ /^-?\d+$/)
         ) {
         $new_level = $getopts_new_level;
      }
      elsif (
          (!defined $new_level) ||
          ($new_level !~ /^-?\d+$/)
         ) {
         $new_level = 0;
         $self->logLocal("No or invalid parameter to setDebugLevel(), " .
                         "setting debug level to 0", 3);
      }

      if ($new_level < 0) {
         $new_level = -1;
      }

      $self->{debug_level} = $new_level;
      $self->logLocal("Set debug level to " . $self->getDebugLevel(), 2);
   }


=item $level = $obj_instance->getDebugLevel();

This function returns the current debug level.  If the current debug
level is not defined, this function returns undefined.

=cut


   sub getDebugLevel() {
      my $self = shift;
      return $self->{debug_level};
   }


=item $obj_instance->setLogFile($log_file);

This function sets the log file name for the C<logLocal()> function.
B<The programmer should call this function before invoking C<setDebugLevel()>>
if the default log file is not to be used.  The function takes one parameter,
C<$log_file>, which defines the new log file name.  If a log file is already
open, it is closed.  The old log file is not truncated or deleted.
Future calls to C<logLocal()> or C<bail()> will log to C<$log_file> if it
is successfully opened.  If the new log file is not successfully opened, 
the function will try to open the default log file, F<program_name.log>.
If that file cannot be opened, F</tmp/program_name.$process_id.log> will
be used.  If no log file argument is passed, the function will try to open
the default log file.  This function is C<GetOptions()> aware; it will accept
two parameters, using the second one as the log file and ignoring the first if
and only if two parameters are passed.  Any other usage specifies the first
parameter as the log file name.

=cut


   sub setLogFile($;$) {
      my $self = shift;
      my $old_log_file = defined $self->{static_log_file} ? 
         $self->{static_log_file} : undef;
      $self->{static_log_file} = shift;
      if (scalar(@_) == 1) {
         $self->{static_log_file} = shift;
      }

      # only consider a new log file that is definable as a file
      if ((defined ($self->{static_log_file})) &&
          ($self->{static_log_file} !~ /^\s*$/)) {
         # delete an old log file entry added by "setLogFile"
         for (my $idx = 0;
              ($idx <= $#{$self->{log_files}}) && defined($old_log_file);
              $idx++) {
            if ($self->{log_files}[$idx] eq $old_log_file) {
               splice @{$self->{log_files}}, $idx, 1;
               $old_log_file = undef;
            }
         }
         unshift @{$self->{log_files}}, $self->{static_log_file};

         # initialize the log file variables and file spaces
         $self->{msg_file_used} = 0;
         $self->{error_file_used} = 0;
         $self->cleanLogFILES();
      }
   }


=item $log_file_name = $obj_instance->getLogFile();

This function returns the name of the log file to be used for printing
log messages.  If no log file is available, this function returns undefined.

=cut


   sub getLogFile() {
      my $self = shift;
      my $return_val = undef;
      if (
          (scalar(@{$self->{log_files}}) != 0) &&
          (defined($self->{log_files}[0]))
         ) {
         $return_val = $self->{log_files}[0];
      }
      return $return_val;
   }


=item $error_file_name = $obj_instance->getErrorFile();

This function returns the name of the error file to be used for printing
error messages.  The error file is derived from the log file; a F<.log>
extension is replaced by a F<.error> extension.  If there is no F<.log>
extension, then F<.error> is appended to the log file name.  If no
log files are defined, this function returns undefined.

=cut


   sub getErrorFile() {
      my $self = shift;
      my $return_val = $self->getLogFile();
      if (defined $return_val) {
         $return_val =~ s/\.log$//g;
         $return_val .= '.error';
      }
      return $return_val;
   }


   # the following private functions are used for logging


   # push items onto the debug level stack
   sub debugPush() {
      my $self = shift;
      if (defined ($self->{debug_level})) {
         push @{$self->{debug_store}}, $self->{debug_level};
      }
      else {
         push @{$self->{debug_store}}, "undef";
      }
      $self->{debug_level} = undef;
   } 


   # pop items from the debug level stack
   sub debugPop() {
      my $self = shift;
      $self->{debug_level} = pop @{$self->{debug_store}};
      if (
          (!defined ($self->{debug_level})) || 
          ($self->{debug_level} eq "undef")
         ) {
         $self->{debug_level} = undef;
      }
   }


   # remove log files
   sub removeLogERROR() {
     
      my $self = shift;
      $self->debugPush();
      if (
          (defined $self->getErrorFile()) &&
          (isWritableFile($self->getErrorFile()))
         ) {
	    unlink $self->getErrorFile() or
            $self->logLocal("Unable to remove error file " .
                             $self->getErrorFile(), 3);
      }
      $self->debugPop();
   }


   sub removeLogMSG() {
      my $self = shift;
      $self->debugPush();
     
      if (
          (defined $self->getLogFile()) &&
          (isWritableFile($self->getLogFile()))
         ) {
            unlink $self->getLogFile() or 
            $self->logLocal("Unable to remove error file " .
                             $self->getLogFile(), 3);
      }
      $self->debugPop();
   }


   # invalidate log files
   sub invalidateLogFILES() {
      my $self = shift;
      $self->debugPush();
      if (defined $self->getLogFile()) {
         $self->logLocal("Invalidating " . $self->getLogFile(), 2); 
         shift @{$self->{log_files}};
         $self->{msg_append_flag} = $self->{error_append_flag} =
            $self->{log_append_setting};
         $self->{msg_file_used} = $self->{error_file_used} = 0;
         $self->cleanLogFILES();
      }
      $self->debugPop();
   }


   # clean previous log files
   sub cleanLogFILES() {
      my $self = shift;
      if ($self->{log_append_setting} == 0) {
         if ($self->{msg_file_used} == 0) {
            $self->removeLogMSG();
         }
         if ($self->{error_file_used} == 0) {
            $self->removeLogERROR();
         }
      }
   }


   # close log files
   sub closeLogERROR() {
      my $self = shift;
      my $return_code = 1; # need to return true for success, false for fail

      $self->debugPush();
      if (!close(ERRLOG) && (defined $self->getErrorFile())) {
         $self->logLocal("Cannot close " . $self->getErrorFile(), 3);
         $return_code = 0;
      }
      else {
         $return_code = 1;
      }
      $self->{error_file_open_flag} = 0;
      $self->debugPop();
      return $return_code;
   }   


   sub closeLogMSG() {
      my $self = shift;
      my $return_code = 1; # need to return true for success, false for fail

      $self->debugPush();
      if (!close(MSGLOG) && (defined $self->getLogFile())) {
         $self->logLocal("Cannot close " . $self->getLogFile(), 3);
         $return_code = 0;
      }
      else {
         $return_code = 1;
      }
      $self->{msg_file_open_flag} = 0;
      $self->debugPop();
      return $return_code;
   }   


   # open log files
   sub openLogERROR() {
      my $self = shift;
      my $return_code = 1; # need to return true for success, false for fail

      $self->debugPush();
      if ((defined $self->getErrorFile()) &&
          ($self->{error_file_open_flag} == 0)) {
         my $fileop;
         $self->{error_file_open_flag} = 1;
         if ($self->{error_append_flag} == 0) {
            $fileop = '>';
            $self->{error_append_flag} = 1;
         }
         else {
            $fileop = '>>';
         }
         if (open(ERRLOG, $fileop . $self->getErrorFile())) {
            autoflush ERRLOG 1;
         }
         else {
            $self->logLocal("Cannot open " . $self->getErrorFile() .
               " for logging", 4);
            $self->{error_file_open_flag} = 0;
         }
      }
      $return_code = $self->{error_file_open_flag};
      $self->debugPop();

      # this is 1 if the file stream is open, 0 if not
      return $return_code;
   }


   sub openLogMSG() {
      my $self = shift;
      my $return_code = 1; # need to return true for success, false for fail

      $self->debugPush();
      if ((defined $self->getLogFile()) && ($self->{msg_file_open_flag} == 0)){
         my $fileop;
         $self->{msg_file_open_flag} = 1;
         if ($self->{msg_append_flag} == 0) {
            $fileop = '>';
            $self->{msg_append_flag} = 1;
         }
         else {
            $fileop = '>>';
         }

         if (open(MSGLOG, $fileop . $self->getLogFile())) {
            autoflush MSGLOG 1;
         }
         else {
            $self->logLocal("Cannot open " . $self->getLogFile() . 
                            " for logging", 4);
            $self->{msg_file_open_flag} = 0;
         }
      }
      $return_code = $self->{msg_file_open_flag};
      $self->debugPop();

      # this is 1 if the file stream is open, 0 if not
      return $return_code;
   }


=item $obj_instance->logAppend($log_append_flag);

The C<logAppend()> function takes either C<0> or C<1> as a flag to 
disable or enable log file appending.  By default, log files are 
truncated at the start of program execution or logging.  Error files are
controlled by this variable as well.  Invalid or undefined calls are ignored.
Calling this function with a C<0> argument after the log files have started
to be written may cause them to be truncated undesirably.  This function is
C<GetOptions()> compliant; if 2 and only 2 variables are passed, the second
option is treated as C<$log_append_flag>.

=cut


   sub logAppend($;$) {
      my $self = shift;
      my $log_append_flag = shift;
      if (defined $_[0]) {
         $log_append_flag = shift;
      }
      if (
          (defined ($log_append_flag)) &&
          (($log_append_flag eq "0") ||
           ($log_append_flag eq "1"))
         ) {
         $self->{log_append_setting} = $self->{msg_append_flag} = 
            $self->{error_append_flag} = $log_append_flag;
      }
   } 


=item $obj_instance->logLocal($log_message, $log_level);

The C<logLocal()> function takes two arguments.  The C<$log_message>
argument specifies the message to be written to the log file.  The
C<$log_level> argument specifies the level at which C<$log_message> is printed.
The active level of logging is set via the C<setDebugLevel()> function.
Only messages at C<$log_level> less than or equal to the active debug
level are logged.  The default debug level is undefined.  Note, a trailing
new line, if it exists, is stripped from the log message.

=cut


   sub logLocal($$) {
      my $self = shift;
      my $log_message = shift;
      my $log_level = shift;

      if ((!defined $log_level) || ($log_level =~ /\D/)) {
         $log_level = 1;
      }

      if (defined $log_message) {
         chomp $log_message; # strip end new line, if it exists

         $log_message = getLogfileDate() . $log_message;
         push @{$self->{debug_queue}}, [ $log_message, $log_level ];

         if ((defined ($self->getDebugLevel())) &&
                      ($self->getDebugLevel() > -1)) {
            while (
                   (defined(my $log_record = $self->{debug_queue}[0])) &&
                   (defined($self->getLogFile()))
                  ) {
               ($log_message, $log_level) = @{$log_record};
               if (
                   (
                    ($log_level <= $self->getDebugLevel()) && # debug level
                    ($self->openLogMSG())                  && # check log file
                    (print MSGLOG "$log_message\n")        && # print message
                    ($self->closeLogMSG())                 && # close log file
                    ($self->{msg_file_used} = 1)              # log file used
                   ) ||
                   (
                    ($log_level > $self->getDebugLevel())     # bad dbg level
                   ) 
                  ) {
                  # log message is successfully processed, so shift it off
                  shift @{$self->{debug_queue}};
               }
               else {
                  $self->debugPush();
                  $self->logLocal("Cannot log message \'$log_message\' to " .
                     $self->getLogFile() . " = " .  $!, 9);
                  $self->invalidateLogFILES(); 
                  $self->debugPop();
               }
            }
         }
      }
      else {
         $self->logLocal("logLocal() called without any parameters!",3);
      }

      while ($#{$self->{debug_queue}} >= $self->{max_debug_queue_size}) {
         # expire old entries; this needs to happen if $self->{debug_level}
         # is undefined or there is no writable log file, otherwise the
         # queue could exhaust RAM.
         shift @{$self->{debug_queue}};
      }
   }
   

=item $obj_instance->logError($log_message,$flag);

The C<logError()> function takes two arguments, the second one being optional. 
The C<$log_message> argument specifies the message to be written to the error 
file. If the C<$flag> argument is defined and is non-zero, the C<$log_message>
is also written to STDERR. The C<$log_message> is also passed to C<logLocal>.
A message passed via logError() will always get logged to the log file 
regardles of the debug level.  Note, a trailing new line, if it exists, is 
stripped from the log message.

=cut


   sub logError($;$) {
      
      my $self = shift;
      my $log_message = shift;
      my $flag = shift; 
      if (defined $log_message) {
         chomp $log_message;  # strip end new line, if it exists
         $self->logLocal($log_message, 0);

         #printing error message to STDERR if flag is non zero.
         if((defined($flag)) && ($flag ne '0')) {
            print STDERR "$log_message\n";
         }        

         $log_message = getLogfileDate() . $log_message;
         push(@{$self->{error_queue}}, $log_message);
        
         while (
                (defined(my $log_message = $self->{error_queue}[0])) &&
                (defined($self->getErrorFile()))
               ) {
               
            if (
                ($self->openLogERROR()) &&     
                (print ERRLOG "$log_message\n") &&
                ($self->closeLogERROR()) &&
                ($self->{error_file_used} = 1) # that is an '='
               ) {
	       shift @{$self->{error_queue}};
            }
            else {
               $self->debugPush();
               $self->logLocal("Cannot log message \'$log_message\' to " .
                  $self->getErrorFile() . " = $!", 6);
               $self->invalidateLogFILES(); 
               $self->debugPop();
            }
         }
      }
      else {
         $self->logLocal("logError() called without any parameters!",3);
      }

      while ($#{$self->{error_queue}} >= $self->{max_debug_queue_size}) {
         # expire old entries; this needs to happen if $self->{debug_level}
         # is undefined or there is no writable log file, otherwise the
         # queue could exhaust RAM.
         shift @{$self->{error_queue}};
      }
   }
 

=item $obj_instance->bail($log_message);

The C<bail()> function takes a single required argument.  The C<$log_message>
argument specifies the message to be passed to C<logLocal()> and displayed
to the screen in using the C<warn> function.  All messages passed to C<bail()>
are logged regardless of the debug level.  The C<bail()> function
calls C<exit(1)> to terminate the program.  Optionally, a second positive
integer argument can be passed as the exit code to use.  Note, a trailing
new line, if it exists, is stripped from the end of the line.

=cut


   sub bail($;$) {
      my $self = shift;
      my $log_message = shift;
      my $exit_code = shift;

      if (
          (!defined $exit_code) ||
          ($exit_code !~ /^\d+$/)
         ) {
         $exit_code = 1;
      }
      if (defined $log_message) {
         chomp $log_message;  # strip end new line, if it exists

         $self->logError($log_message);
         print STDERR $log_message, "\n";
      }

      exit $exit_code;
   }


# Functional Class : modified methods

=item $getopts_error_code = $obj_instance->TIGR_GetOptions(@getopts_arguments);

This function extends C<Getopt::Long::GetOptions()>.  It may be used
as C<GetOptions()> is used.  Extended functionality eliminates the need
to C<eval {}> the block of code containing the function.  Further, TIGR
standard options, such as C<-help>, are defined implicitly.  Using this
function promotes proper module behavior.  Log and error files from 
previous runs are removed if the log file append option, C<-appendlog>,
is not set to 1.

The following options are defined by this function:

=over

=item -appendlog APPEND_FLAG

Passing '1' to this argument turns on log file appending.

=item -debug DEBUG_LEVEL

Set debugging to DEBUG_LEVEL.

=item -logfile LOG_FILE_NAME

Set the default TIGR Foundation log file to LOG_FILE_NAME.

=item -version, -V

Print version information and exit.

=item -help, -h

Print help information and exit.

=item -depend

Print dependency information and exit.

=back

Regular C<GetOptions()> may still be used, however the C<TIGR_GetOptions()>
function eliminates some of the confusing issues with setting log files
and debug levels.  B<The options defined by C<TIGR_GetOptions()> cannot be
overridden or recorded>.  To get the log file and debug level after parsing
the command line, use C<getLogFile()> and C<getDebugLevel()>.  C<GetOptions()>
default variables, ie. those of the form C<$opt_I<optionname>>, are not
supported.  This function will return 1 on success.

=cut


   sub TIGR_GetOptions(@) {
      my $self = shift;
      my @user_options = @_;

      my $appendlog_var = undef;
      my $logfile_var = undef;
      my $debug_var = undef;
      my $version_var = undef;
      my $help_var = undef;
      my $depend_var = undef;

      # these foundation options support the defaults
      my @foundation_options = (
         "appendlog=i" => \$appendlog_var,
         "logfile=s" => \$logfile_var,
         "debug=i" => \$debug_var,
         "version|V" => \$version_var,
         "help|h" => \$help_var,
         "depend" => \$depend_var
      );

      Getopt::Long::Configure('no_ignore_case');
      my $getopt_code = eval 'GetOptions (@user_options, @foundation_options)';

      if ((defined $help_var) && ($help_var =~ /^(.*)$/))  {
         $self->printHelpInfoAndExit();
      }

      if ((defined $version_var) && ($version_var =~ /^(.*)$/)) {
         $self->printVersionInfoAndExit();
      }

      if ((defined $depend_var) && ($depend_var =~ /^(.*)$/)) {
         $self->printDependInfoAndExit();
      }

      if ((defined $appendlog_var) && ($appendlog_var =~ /^(.*)$/)) {
	 $appendlog_var = $1;
         $self->logAppend($appendlog_var);
      }

      if ((defined $logfile_var) && ($logfile_var =~ /^(.*)$/))  {
	 $logfile_var = $1;
         $self->setLogFile($logfile_var);
      }

      if ((defined $debug_var) && ($debug_var =~ /^(.*)$/)) {
	 $debug_var = $1;
         $self->setDebugLevel($debug_var);
      }

      # remove old log files, if necessary
      for (
           my $file_control_var = 0;
           $file_control_var <= $#{$self->{log_files}};
           $file_control_var++
          ) {
          $self->cleanLogFILES();
          push(@{$self->{log_files}}, shift @{$self->{log_files}});
      }
      return $getopt_code;
   }

   DESTROY {
      my $self = shift;
      $self->{finish_time} = time;
      my $time_difference = $self->{finish_time} - $self->{start_time};
      my $num_days = int($time_difference / 86400); # there are 86400 sec/day
      $time_difference -= $num_days * 86400;
      my $num_hours = int($time_difference / 3600); # there are 3600 sec/hour
      $time_difference -= $num_hours * 3600;
      my $num_min = int($time_difference / 60); # there are 60 sec/hour
      $time_difference -= $num_min * 60;
      my $num_sec = $time_difference; # the left overs are seconds
      my $time_str = sprintf "%03d-%02d:%02d:%02d", $num_days, $num_hours,
         $num_min, $num_sec;
      $self->logLocal("FINISH: " . $self->getProgramInfo('name') .
         ", elapsed ".$time_str ,0);
   }
}

=back

=head1 USAGE

To use this module, load the C<TIGR::Foundation> package
via the C<use> function.  Then, create a new instance of the
object via the C<new()> method, as shown below.  If applicable,
C<START> and C<FINISH> log messages are printed when the object
is created and destroyed, respectively.  It is advisable to 
keep the instance of the object in scope for the whole program
to achieve maximum functionality.

An example script using this module follows:

   use strict;
   use TIGR::Foundation;

   my $tfobject = new TIGR::Foundation;

   MAIN: 
   {
      # The following dependencies are not used in
      # this script, but are provided as an example.

      my @DEPEND = ("/usr/bin/tee", "/sbin/stty");

      # The user defined $VERSION variable is usable by Perl.
      # The auto defined $REVISION variable stores the RCS/CVS revision
      # The user defined $VERSION_STRING reports both.

      my $VERSION = '1.0';
      my $REVISION = (qw$Revision: 1.1.1.1 $)[-1];
      my $VERSION_STRING = "$VERSION (Build $REVISION)";

      my $HELP_INFO = "This is my help\n";

      # All of the necessary information must be passed
      # to the foundation object instance, as below.

      $tfobject->addDependInfo(@DEPEND);
      $tfobject->setVersionInfo($VERSION_STRING);
      $tfobject->setHelpInfo($HELP_INFO);

      my $input_file;
      my $output_file;

      $tfobject->TIGR_GetOptions("input=s" => \$input_file,
                                 "output=s" => \$output_file);

      # GetOptions(), and subsequently TIGR_GetOptions(), leaves
      # the variables unchanged if no corresponding command line
      # arguments are parsed.  The passed variables are checked below.

      if (defined $input_file) {

         # The log message is written only if debugging is turned on.
         # By default, debugging is off.  To turn on debugging, use the
         # '-debug DEBUG_LEVEL' option on the command line.
         # In this example, '-debug 1' would set debugging to level 1
         # and report these log messages.

         $tfobject->logLocal("My input file is $input_file", 1);
      }

      print "Hello world", "\n";

      # This case is similar to the previous one above...
      if (defined $output_file) {
         $tfobject->logLocal("My output file is $output_file.", 1);
      }
   }

=cut

1; 




