<embed src="../specialization_project.pdf"
       style="position: fixed; 
              left: 50%;
              transform: translate(-50%, 0%);
              width: 100vw;"/>

<script>
    header = document.getElementsByTagName("header")[0];
    pdf = document.getElementsByTagName("embed")[0];
    offset = header.offsetHeight;
    pdfHeight = window.innerHeight - offset;
    pdf.style.top = offset + "px";
    pdf.style.height = (pdfHeight ) + "px";

</script>