$(document).ready(function() { 
    $('#gene-search').submit(function (event) {
        window.location.href = "/gene#" + $('#search-input').val();
        return false; // prevents form submission from overriding redirect
    });

    $.ajax({
        dataType: "json",
        url: "/data/rn7.autocomplete.json"
    }).done(function (data) {
        $("#search-input").autocomplete({
            // source: data,
            source: function (request, response) {
                var results = $.ui.autocomplete.filter(data, request.term);
                response(results.slice(0, 10));
            },
            select: function(event, ui) {
                if (ui.item) {
                    $(this).val(ui.item.value);
                    $(this.form).submit();
                }
            }
        });
    });
});
