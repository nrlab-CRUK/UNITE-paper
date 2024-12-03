#!/usr/bin/perl
use strict;
use warnings;
use File::Find;

# Define the subroutine to rename the files
sub rename_files {
    my $file = $_;
    
    # Only proceed if the filename contains 'xgboost_'
    if ($file =~ /xgboost_/) {
        my $new_name = $file;
        $new_name =~ s/xgboost_/lr_/g;  # Replace 'xgboost_' with 'lr_'
        
        # Ensure we're renaming a file and not a directory
        if (-f $file) {
            my $dir = $File::Find::name;
            $dir =~ s|$file$||;  # Extract the directory part from the full path
            
            # Rename the file
            rename($file, "$dir$new_name") or warn "Could not rename '$file' to '$new_name': $!\n";
            print "Renamed '$file' to '$new_name'\n";
        }
    }
}

# Start the recursive search from the current directory
find(\&rename_files, '.');
