# following documentation from Microsoft
# https://docs.microsoft.com/en-us/azure/virtual-machines/linux/mount-azure-file-storage-on-linux-using-smb


export STORAGEACCT="h2ab3data"
STORAGEKEY=$(az storage account keys list \
    --resource-group "H2AB3" \
    --account-name $STORAGEACCT \
    --query "[0].value" | tr -d '"')


az storage file upload-batch \
    --account-name $STORAGEACCT \
    --account-key $STORAGEKEY \
    --source "RNA-Seq" \
    --destination-path "." \
    --destination "https://h2ab3data.file.core.windows.net/raw"


# creating a mountpoint on the DSVM and mounting azure files shares
STORAGEKEY="SoNLLyJ0ci0uFfJZgX+qJ7qAPO82bPliyOv2Pg79UzWm41Zn8sTxCO8uuxl6psVgPjtCYE0KLOaNw3IPGuOvfQ=="
DEST="h2ab3data.file.core.windows.net/raw"
sudo mount -t cifs //$DEST /mnt/raw -o vers=3.0,username=$STORAGEACCT,password=$STORAGEKEY,dir_mode=0777,file_mode=0777,serverino

DEST="h2ab3data.file.core.windows.net/processed"
sudo mount -t cifs //$DEST /mnt/processed -o vers=3.0,username=$STORAGEACCT,password=$STORAGEKEY,dir_mode=0777,file_mode=0777,serverino

# /etc/fstab entries to make it persistent
//h2ab3data.file.core.windows.net/processed /mnt/processed cifs vers=3.0,username=h2ab3data,password=SoNLLyJ0ci0uFfJZgX+qJ7qAPO82bPliyOv2Pg79UzWm41Zn8sTxCO8uuxl6psVgPjtCYE0KLOaNw3IPGuOvfQ==,dir_mode=0777,file_mode=0777
//h2ab3data.file.core.windows.net/raw /mnt/raw cifs vers=3.0,username=h2ab3data,password=SoNLLyJ0ci0uFfJZgX+qJ7qAPO82bPliyOv2Pg79UzWm41Zn8sTxCO8uuxl6psVgPjtCYE0KLOaNw3IPGuOvfQ==,dir_mode=0777,file_mode=0777
